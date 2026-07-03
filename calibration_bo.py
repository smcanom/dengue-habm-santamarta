#!/usr/bin/env python3
"""
Surrogate / Bayesian-optimization calibration of the Santa Marta dengue HABM.

Calibrates four parameters jointly:
    I0       initial infected humans  (count.of.infected.humans)
    z        human-interaction parameter
    beta_mh  mosquito->human transmission probability
    w        REPORTING DELAY in weeks  (post-hoc shift of the simulated output)

Why a reporting delay (not a longer run): a case with *onset* in simulated week j
is *reported* in week j+w. We therefore run the normal 52-week (365-tick) simulation
-- fully inside the one year of temperature data -- and shift the simulated weekly
onset series later by w when scoring it against the observed *reported* series. This
needs no horizon extension and no temperature wrap-around.

Search uses scikit-optimize (GP surrogate + Expected Improvement), evaluating points
in parallel batches. Each simulation runs in an ISOLATED working directory (the model
reads/writes hard-coded CSV names in the CWD), so jobs never clobber one another.

Objective (equal-weight composite, computed on the w-shifted overlap):
    level  = RMSE / std(obs)                 (magnitude agreement)
    timing = |argmax(sim) - argmax(obs)| / L  (residual peak-timing error)
    shape  = (1 - Pearson(zscore(smooth(sim)), zscore(smooth(obs)))) / 2

Usage (inside the uv env):
    uv run python calibration_bo.py --smoke            # fast plumbing test
    uv run python calibration_bo.py --n-calls 40 --replicates 4 --max-workers 14
"""
from __future__ import annotations
import argparse
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

# Reuse vetted helpers from the existing I0 calibration.
from calibration_I0 import rewrite_props, read_results_csv, daily_to_weekly, read_observed

# Input files the model reads (by hard-coded relative name) from its CWD.
# NB: do NOT include output files (results.csv, weekly_neighborhood_cases.csv,
# neighborhood_cells.csv, summary_neighborhood_sizes.csv, hi_activation_map.csv) --
# symlinking those would make the model write back into the repo root.
INPUT_FILES = [
    "santa_marta_grid_relative.csv",
    "SantaMartaChange.csv",
    "temperatura-precipitacion-bello.csv",
    "hi_by_neighborhood_2023.csv",
    "locations.csv",
]

REPAST_LIB_DEFAULT = "/home/scano/sfw/repast_hpc-2.3.1/lib"


# ----------------------------------------------------------------------------
# Simulation execution (isolated, cached)
# ----------------------------------------------------------------------------

def _sim_env(repast_lib: str) -> Dict[str, str]:
    env = dict(os.environ)
    prev = env.get("LD_LIBRARY_PATH", "")
    env["LD_LIBRARY_PATH"] = repast_lib + (":" + prev if prev else "")
    env["OMPI_MCA_btl_vader_single_copy_mechanism"] = "none"
    env["OMP_NUM_THREADS"] = "1"  # avoid thread oversubscription under the pool
    return env


def _cache_key(I0: int, z: float, beta_mh: float, n_humans: int, stop_at: int, seed: int) -> str:
    return f"I0-{I0}_z-{z:.4f}_bmh-{beta_mh:.4f}_N-{n_humans}_T-{stop_at}_seed-{seed:06d}"


def run_one_sim(args, I0: int, z: float, beta_mh: float, seed: int) -> np.ndarray:
    """Run ONE replicate in an isolated dir; return the weekly onset series.

    Results are cached on disk by parameter+seed key, so resumes and duplicate
    proposals never recompute.
    """
    key = _cache_key(I0, z, beta_mh, args.n_humans, args.stop_at, seed)
    cache = args.outdir / "sims" / f"{key}.npy"
    if cache.exists():
        return np.load(cache)

    rundir = args.rundir / key
    rundir.mkdir(parents=True, exist_ok=True)
    # Symlink inputs into the isolated working dir.
    for fname in INPUT_FILES:
        src = (args.workdir / fname).resolve()
        dst = rundir / fname
        if src.exists() and not dst.exists():
            dst.symlink_to(src)

    # Per-job props: copy the nominal baseline, then override the calibrated values.
    job_props = rundir / "model.props"
    shutil.copy(args.baseline_props, job_props)
    rewrite_props(job_props, {
        "count.of.humans": str(args.n_humans),
        "count.of.infected.humans": str(int(I0)),
        "stop.at": str(args.stop_at),
        "z": f"{z:.6f}",
        # Corrected model: mosquito<-human infection uses beta_hm (paper Eq. psi/mosquito_dynamics).
        # Tie beta_hm = beta_mh = beta so a single transmission knob is calibrated.
        "beta_mh": f"{beta_mh:.6f}",
        "beta_hm": f"{beta_mh:.6f}",
        "random.seed": str(seed),
        "global.random.seed": str(seed),
        "base_seed": str(args.base_seed),
        "replicate_id": str(seed - args.base_seed),
        "config_id": "bo",
        "perturb_param": "baseline",
        "perturb_delta": "0.0",
        "obs_csv": str(args.observed.resolve()),
    })

    # --bind-to none / --oversubscribe: each independent `mpirun -n 1` otherwise
    # binds its rank to core 0, so concurrent jobs all pile onto cores 0-1 and
    # crawl (~36x slower). Letting the OS schedule them spreads across all cores.
    cmd = ["mpirun", "--bind-to", "none", "--oversubscribe", "-n", "1",
           str(args.main_exe.resolve()),
           str(args.config_props.resolve()), "model.props"]
    proc = subprocess.run(cmd, cwd=str(rundir), env=_sim_env(args.repast_lib),
                          stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    res = rundir / "results.csv"
    if proc.returncode != 0 or not res.exists():
        err = proc.stderr.decode("utf-8", "ignore")[-500:]
        raise RuntimeError(f"sim failed ({key}) rc={proc.returncode}: {err}")

    weekly = daily_to_weekly(read_results_csv(res), start_tick=1)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.save(cache, weekly)
    if not args.keep_rundirs:
        shutil.rmtree(rundir, ignore_errors=True)
    return weekly


# ----------------------------------------------------------------------------
# Composite objective (reporting-delay aware)
# ----------------------------------------------------------------------------

def _smooth_z(x: np.ndarray) -> np.ndarray:
    s = pd.Series(x).rolling(3, center=True, min_periods=1).mean().to_numpy()
    return (s - s.mean()) / (s.std() + 1e-12)


def composite_score(sim_onset: np.ndarray, obs: np.ndarray, w: int, L: int
                    ) -> Tuple[float, Dict[str, float]]:
    """Equal-weight composite on the reporting-delay-shifted overlap.

    Onset week j is reported at week j+w, so we compare observed[w:L] against
    simulated onset[0:L-w].
    """
    m = L - int(w)
    if m < 20:  # too little overlap to score reliably
        return 5.0, {"level": 5.0, "timing": 1.0, "shape": 1.0, "rmse": float("nan")}
    sim = np.asarray(sim_onset, float)[:m]
    o = np.asarray(obs, float)[int(w):int(w) + m]
    n = min(len(sim), len(o))
    sim, o = sim[:n], o[:n]

    rmse = float(np.sqrt(np.mean((sim - o) ** 2)))
    level = rmse / (o.std() + 1e-9)
    timing = abs(int(np.argmax(sim)) - int(np.argmax(o))) / float(n)
    zs, zo = _smooth_z(sim), _smooth_z(o)
    pear = float(np.corrcoef(zs, zo)[0, 1]) if n > 2 else float("nan")
    shape = 1.0 if np.isnan(pear) else (1.0 - pear) / 2.0
    score = (level + timing + shape) / 3.0
    return score, {"level": level, "timing": timing, "shape": shape,
                   "rmse": rmse, "pearson": pear}


def composite_score_bestwindow(sim_onset: np.ndarray, obs_full: np.ndarray, L: int
                               ) -> Tuple[float, Dict[str, float]]:
    """Score an L-week simulated season against the best-matching L-week window
    of the longer (multi-year) observed record.

    The observed record spans 2023-2024, so it contains more than one possible
    one-year window. For each candidate we slide a length-L frame across the
    record and keep the window with the lowest composite score. The chosen
    offset is the phase alignment; it replaces the old reporting-delay search.
    """
    sim = np.asarray(sim_onset, float)
    Ls = min(int(L), len(sim))
    sim = sim[:Ls]
    obs_full = np.asarray(obs_full, float)
    max_off = max(0, len(obs_full) - Ls)
    best_score = 5.0
    best = {"level": 5.0, "timing": 1.0, "shape": 1.0, "rmse": float("nan"),
            "pearson": float("nan"), "offset": 0}
    if Ls < 20:
        return best_score, best
    for off in range(0, max_off + 1):
        o = obs_full[off:off + Ls]
        if len(o) < Ls:
            break
        rmse = float(np.sqrt(np.mean((sim - o) ** 2)))
        level = rmse / (o.std() + 1e-9)
        timing = abs(int(np.argmax(sim)) - int(np.argmax(o))) / float(Ls)
        zs, zo = _smooth_z(sim), _smooth_z(o)
        pear = float(np.corrcoef(zs, zo)[0, 1]) if Ls > 2 else float("nan")
        shape = 1.0 if np.isnan(pear) else (1.0 - pear) / 2.0
        score = (level + timing + shape) / 3.0
        if score < best_score:
            best_score = score
            best = {"level": level, "timing": timing, "shape": shape,
                    "rmse": rmse, "pearson": pear, "offset": off}
    return best_score, best


# ----------------------------------------------------------------------------
# Batch evaluation: dedup sims, run in parallel, assemble per-point scores
# ----------------------------------------------------------------------------

def evaluate_points(args, pts: List[list], observed: np.ndarray, L: int, pool: ThreadPoolExecutor
                    ) -> List[Tuple[float, dict]]:
    seeds = [args.base_seed + r for r in range(args.replicates)]
    # Unique (I0,z,beta,seed) jobs across the batch.
    jobs: Dict[str, Tuple[int, float, float, int]] = {}
    for (I0, z, beta_mh) in pts:
        for s in seeds:
            jobs[_cache_key(int(I0), z, beta_mh, args.n_humans, args.stop_at, s)] = (int(I0), z, beta_mh, s)

    futs = {k: pool.submit(run_one_sim, args, *v) for k, v in jobs.items()}
    results = {k: f.result() for k, f in futs.items()}

    out = []
    for (I0, z, beta_mh) in pts:
        reps = [results[_cache_key(int(I0), z, beta_mh, args.n_humans, args.stop_at, s)] for s in seeds]
        T = min(len(r) for r in reps)
        mean_onset = np.vstack([r[:T] for r in reps]).mean(axis=0)
        score, comp = composite_score_bestwindow(mean_onset, observed, L)
        comp.update({"I0": int(I0), "z": float(z), "beta_mh": float(beta_mh)})
        out.append((score, comp))
    return out


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description="Bayesian-optimization calibration {I0,z,beta_mh,w}.")
    ap.add_argument("--observed", type=Path, default=Path("view/incidence_rates_santa_marta.csv"))
    ap.add_argument("--baseline-props", type=Path, default=Path("props/model_baseline.props"))
    ap.add_argument("--config-props", type=Path, default=Path("props/config.props"))
    ap.add_argument("--main-exe", type=Path, default=Path("./Main.exe"))
    ap.add_argument("--repast-lib", type=str, default=REPAST_LIB_DEFAULT)
    ap.add_argument("--workdir", type=Path, default=Path("."))
    ap.add_argument("--outdir", type=Path, default=Path("calib_bo_runs"))
    ap.add_argument("--rundir", type=Path, default=Path("_bo_workdirs"))
    ap.add_argument("--n-humans", type=int, default=100000)
    ap.add_argument("--stop-at", type=int, default=365)
    ap.add_argument("--weeks", type=int, default=52, help="L: comparison-window length in weeks.")
    ap.add_argument("--warmup-max", type=int, default=12, help="Max reporting delay w (weeks).")
    ap.add_argument("--n-calls", type=int, default=40)
    ap.add_argument("--n-initial", type=int, default=10)
    ap.add_argument("--batch", type=int, default=4, help="BO points proposed per iteration.")
    ap.add_argument("--replicates", type=int, default=4, help="CRN replicates per evaluation.")
    ap.add_argument("--final-replicates", type=int, default=15)
    ap.add_argument("--max-workers", type=int, default=min(14, os.cpu_count() or 4))
    ap.add_argument("--base-seed", type=int, default=12345)
    ap.add_argument("--z-bounds", type=float, nargs=2, default=[0.05, 0.5])
    # Corrected model: R0 crosses 1 near beta~0.20-0.25; broad waves occur just above
    # threshold. Search the near-critical band (beta_mh=beta_hm=beta tied in run_one_sim).
    ap.add_argument("--bmh-bounds", type=float, nargs=2, default=[0.19, 0.32])
    ap.add_argument("--I0-bounds", type=int, nargs=2, default=[2, 8])
    ap.add_argument("--seed", type=int, default=0, help="BO random_state.")
    ap.add_argument("--evaluate-only", type=str, default=None,
                    help='Skip BO; run final-replicates at a fixed "I0,z,beta_mh,w" and write an ensemble.')
    ap.add_argument("--ensemble-out", type=Path, default=None,
                    help="Output CSV for the ensemble (defaults to outdir/final_ensemble_bo.csv).")
    ap.add_argument("--keep-rundirs", action="store_true")
    ap.add_argument("--smoke", action="store_true", help="Fast plumbing test (tiny scale).")
    args = ap.parse_args()

    if args.smoke:
        args.n_humans = 3000
        args.stop_at = 56          # 8 weeks
        args.weeks = 8
        args.warmup_max = 2
        args.n_calls = 4
        args.n_initial = 3
        args.batch = 2
        args.replicates = 2
        args.final_replicates = 3
        args.max_workers = 4

    from skopt import Optimizer
    from skopt.space import Real, Integer

    args.outdir.mkdir(parents=True, exist_ok=True)
    args.rundir.mkdir(parents=True, exist_ok=True)

    weeks_obs, y = read_observed(args.observed)
    # observed is the FULL multi-year record; L is the one-year comparison length.
    observed = np.asarray(y, float)
    L = int(args.weeks)
    print(f"[BO] observed weeks={len(observed)} (full record) -> compare L={L}; "
          f"N={args.n_humans} stop.at={args.stop_at} "
          f"reps={args.replicates} n_calls={args.n_calls} workers={args.max_workers}")

    # ---- Fixed-point mode (e.g., nominal baseline): no BO, just an ensemble. ----
    if args.evaluate_only:
        a = [x.strip() for x in args.evaluate_only.split(",")]
        I0f, zf, bf, wf = int(round(float(a[0]))), float(a[1]), float(a[2]), int(round(float(a[3])))
        out_path = args.ensemble_out or (args.outdir / "final_ensemble_bo.csv")
        seeds = [args.base_seed + 10_000 + r for r in range(args.final_replicates)]
        print(f"[EVAL] fixed point I0={I0f} z={zf} beta_mh={bf} w={wf}, {args.final_replicates} reps")
        with ThreadPoolExecutor(max_workers=args.max_workers) as pool:
            reps = [f.result() for f in
                    [pool.submit(run_one_sim, args, I0f, zf, bf, s) for s in seeds]]
        T = min(len(r) for r in reps)
        W = np.vstack([r[:T] for r in reps])
        score, comp = composite_score(W.mean(axis=0), observed, wf, L)
        print(f"[EVAL] composite={score:.4f} level={comp['level']:.3f} timing={comp['timing']:.3f} "
              f"shape={comp['shape']:.3f} pearson={comp.get('pearson', float('nan')):.3f} rmse={comp['rmse']:.3f}")
        m = min(L - wf, T)
        pd.DataFrame({
            "week_idx": np.arange(m), "obs": observed[wf:wf + m],
            "median": np.median(W, axis=0)[:m], "p15": np.percentile(W, 15, axis=0)[:m],
            "p85": np.percentile(W, 85, axis=0)[:m], "mean": W.mean(axis=0)[:m],
        }).to_csv(out_path, index=False)
        print(f"[EVAL] wrote {out_path}")
        return

    space = [Integer(args.I0_bounds[0], args.I0_bounds[1], name="I0"),
             Real(args.z_bounds[0], args.z_bounds[1], name="z"),
             Real(args.bmh_bounds[0], args.bmh_bounds[1], name="beta_mh")]
    opt = Optimizer(space, base_estimator="GP", acq_func="EI",
                    n_initial_points=args.n_initial, random_state=args.seed)

    trials: List[dict] = []
    t0 = time.time()
    done = 0
    with ThreadPoolExecutor(max_workers=args.max_workers) as pool:
        while done < args.n_calls:
            b = min(args.batch, args.n_calls - done)
            pts = opt.ask(n_points=b)
            evals = evaluate_points(args, pts, observed, L, pool)
            opt.tell(pts, [e[0] for e in evals])
            for (score, comp) in evals:
                comp["score"] = score
                trials.append(comp)
                done += 1
                print(f"  [{done:3d}/{args.n_calls}] score={score:.4f} "
                      f"I0={comp['I0']} z={comp['z']:.3f} bmh={comp['beta_mh']:.3f} "
                      f"offset={comp['offset']} (level={comp['level']:.3f} timing={comp['timing']:.3f} "
                      f"shape={comp['shape']:.3f})")

    df = pd.DataFrame(trials).sort_values("score").reset_index(drop=True)
    df.to_csv(args.outdir / "bo_trials.csv", index=False)
    best = df.iloc[0]
    print(f"\n[BO] best score={best['score']:.4f} | I0={int(best['I0'])} z={best['z']:.4f} "
          f"beta_mh={best['beta_mh']:.4f} window_offset={int(best['offset'])}  ({time.time()-t0:.0f}s)")

    # Final ensemble at the best point (more replicates) for median + 70% band.
    print(f"[BO] final ensemble: {args.final_replicates} replicates at best point...")
    seeds = [args.base_seed + 10_000 + r for r in range(args.final_replicates)]
    with ThreadPoolExecutor(max_workers=args.max_workers) as pool:
        futs = [pool.submit(run_one_sim, args, int(best["I0"]), float(best["z"]),
                            float(best["beta_mh"]), s) for s in seeds]
        reps = [f.result() for f in futs]
    T = min(len(r) for r in reps)
    W = np.vstack([r[:T] for r in reps])
    # Re-select the best-matching window for the final ensemble mean.
    _, fcomp = composite_score_bestwindow(W.mean(axis=0), observed, L)
    off = int(fcomp["offset"])
    m = min(L, T, len(observed) - off)
    idx = np.arange(m)
    out = pd.DataFrame({
        "week_idx": idx,
        "obs": observed[off:off + m],
        "median": np.median(W, axis=0)[:m],
        "p15": np.percentile(W, 15, axis=0)[:m],
        "p85": np.percentile(W, 85, axis=0)[:m],
        "mean": W.mean(axis=0)[:m],
    })
    out.to_csv(args.outdir / "final_ensemble_bo.csv", index=False)
    print(f"[BO] best window offset={off}  Pearson={fcomp.get('pearson', float('nan')):.3f}")
    print(f"[BO] wrote {args.outdir/'bo_trials.csv'} and {args.outdir/'final_ensemble_bo.csv'}")


if __name__ == "__main__":
    sys.exit(main())
