#!/usr/bin/env python3
"""
Corrected 7-parameter one-at-a-time (OAT) sensitivity analysis for the
paper-aligned dengue HABM.

Why this exists
---------------
The published SA perturbed sigma_M, sigma_H, z, beta_mh -- but in the original
code sigma_M/sigma_H were hard-coded constants and z was applied to the wrong
(mosquito) population, so three of those four knobs were inert. After wiring the
code to the paper's Eqs. (bites), (human_infection), (psi), (mosquito_dynamics),
all transmission parameters are live and this SA re-measures their influence.

Operating point
---------------
The model only matches the observed season *on* its bifurcation (beta~0.28),
where a +10% perturbation flips ~1% -> ~39% attack: SA there measures distance
to the threshold, not parameter influence. So this SA runs at a robustly
super-critical nominal (default beta=0.35, z=0.118) where the epidemic reliably
occurs and perturbations scale it smoothly. Report as "sensitivity in the
epidemic regime", distinct from the (bifurcation) fitting point.

Method
------
One-sided OAT: baseline + each parameter individually increased by +10%, run at
N CRN replicates (seed = base_seed + replicate_id, shared across configs). The
perturbed VALUE is written directly into props with perturb_param=baseline so
Main.cpp does not double-apply its own perturbation. NSC = proportional change
in the mean 70% prediction-band width of weekly incidence (paper's metric); peak
height/timing changes are also reported.
"""
import argparse
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List

import numpy as np

from calibration_I0 import rewrite_props, read_results_csv, daily_to_weekly

INPUT_FILES = [
    "santa_marta_grid_relative.csv",
    "SantaMartaChange.csv",
    "temperatura-precipitacion-bello.csv",
    "hi_by_neighborhood_2023.csv",
    "locations.csv",
]
REPAST_LIB_DEFAULT = "/home/scano/sfw/repast_hpc-2.3.1/lib"

# Robustly super-critical nominal for the corrected model (see module docstring).
NOMINAL = {
    "beta_mh": 0.35,
    "beta_hm": 0.35,
    "z": 0.118264,
    "sigma_M": 0.3,
    "sigma_H": 3.0,
    "r": 0.6,
    "C": 30.0,
}
PERTURB_PARAMS = ["beta_mh", "beta_hm", "z", "sigma_M", "sigma_H", "r", "C"]
PERTURB_DELTA = 0.10


def _sim_env(repast_lib: str) -> Dict[str, str]:
    env = dict(os.environ)
    prev = env.get("LD_LIBRARY_PATH", "")
    env["LD_LIBRARY_PATH"] = repast_lib + (":" + prev if prev else "")
    env["OMPI_MCA_btl_vader_single_copy_mechanism"] = "none"
    env["OMP_NUM_THREADS"] = "1"
    return env


def run_one(args, cfg_id: str, params: Dict[str, float], seed: int) -> np.ndarray:
    """Run one (config, seed) in an isolated working dir; return weekly incidence."""
    rundir = Path(args.workdir) / f"{cfg_id}_seed{seed}"
    if rundir.exists():
        shutil.rmtree(rundir, ignore_errors=True)
    rundir.mkdir(parents=True)
    for name in INPUT_FILES:
        src = (Path.cwd() / name).resolve()
        if src.exists():
            (rundir / name).symlink_to(src)

    job_props = rundir / "model.props"
    shutil.copy(args.baseline_props, job_props)
    updates = {
        "count.of.humans": str(args.n_humans),
        "count.of.infected.humans": str(args.i0),
        "stop.at": str(args.stop_at),
        "random.seed": str(seed),
        "global.random.seed": str(seed),
        "base_seed": str(args.base_seed),
        "replicate_id": str(seed - args.base_seed),
        # write perturbed values directly; keep Main.cpp's own perturbation off
        "perturb_param": "baseline",
        "perturb_delta": "0.0",
        "config_id": cfg_id,
        "obs_csv": str(Path(args.observed).resolve()),
    }
    for k, v in params.items():
        updates[k] = f"{v:.6f}"
    rewrite_props(job_props, updates)

    cmd = ["mpirun", "--bind-to", "none", "--oversubscribe", "-n", "1",
           str(Path(args.main_exe).resolve()),
           str(Path(args.config_props).resolve()), "model.props"]
    proc = subprocess.run(cmd, cwd=str(rundir), env=_sim_env(args.repast_lib),
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (rundir / "run.log").write_bytes(
        b"=== STDOUT ===\n" + proc.stdout + b"\n=== STDERR ===\n" + proc.stderr)
    res = rundir / "results.csv"
    if proc.returncode != 0 or not res.exists():
        err = proc.stderr.decode("utf-8", "ignore")[-400:]
        # leave rundir in place for inspection on failure
        raise RuntimeError(f"sim failed ({cfg_id} seed {seed}) rc={proc.returncode} "
                           f"results.csv={res.exists()}: {err}")
    weekly = daily_to_weekly(read_results_csv(res), start_tick=1)
    if not args.keep_rundirs:
        shutil.rmtree(rundir, ignore_errors=True)
    return weekly


def band_width(M: np.ndarray) -> float:
    """Mean 70% prediction-band width (p85 - p15) across weeks."""
    p15 = np.percentile(M, 15, axis=0)
    p85 = np.percentile(M, 85, axis=0)
    return float(np.mean(p85 - p15))


def main():
    ap = argparse.ArgumentParser(description="Corrected 7-parameter OAT sensitivity analysis.")
    ap.add_argument("--replicates", type=int, default=5)
    ap.add_argument("--base-seed", type=int, default=12345)
    ap.add_argument("--n-humans", type=int, default=100000)
    ap.add_argument("--i0", type=int, default=8)
    ap.add_argument("--stop-at", type=int, default=365)
    ap.add_argument("--max-workers", type=int, default=min(14, os.cpu_count() or 4))
    ap.add_argument("--baseline-props", default="props/model_baseline.props")
    ap.add_argument("--config-props", default="props/config.props")
    ap.add_argument("--main-exe", default="./Main.exe")
    ap.add_argument("--observed", default="view/incidence_rates_santa_marta.csv")
    ap.add_argument("--repast-lib", default=REPAST_LIB_DEFAULT)
    ap.add_argument("--workdir", default="_sa_workdirs")
    ap.add_argument("--outdir", default="sa_corrected_runs")
    ap.add_argument("--keep-rundirs", action="store_true")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    seeds = [args.base_seed + r for r in range(args.replicates)]

    # Build configs: baseline + one +10% perturbation each.
    configs = [("baseline", dict(NOMINAL))]
    for p in PERTURB_PARAMS:
        pert = dict(NOMINAL)
        pert[p] = NOMINAL[p] * (1.0 + PERTURB_DELTA)
        configs.append((f"{p}_p10", pert))

    print(f"[SA] nominal={NOMINAL}")
    print(f"[SA] {len(configs)} configs x {args.replicates} reps "
          f"= {len(configs)*args.replicates} runs, N={args.n_humans} stop.at={args.stop_at} "
          f"workers={args.max_workers}")

    # Launch all (config, seed) jobs concurrently.
    jobs = [(cid, params, s) for (cid, params) in configs for s in seeds]
    results: Dict[tuple, np.ndarray] = {}
    with ThreadPoolExecutor(max_workers=args.max_workers) as ex:
        futs = {ex.submit(run_one, args, cid, params, s): (cid, s) for (cid, params, s) in jobs}
        for fut in futs:
            cid, s = futs[fut]
            try:
                results[(cid, s)] = fut.result()
            except Exception as e:
                print(f"  [FAIL] {cid} seed {s}: {e}", file=sys.stderr)

    # Aggregate per config (align weekly series to min length).
    rows = []
    base_bw = None
    base_peak = None
    per_config_weekly = {}
    for (cid, _params) in configs:
        series = [results[(cid, s)] for s in seeds if (cid, s) in results]
        if not series:
            print(f"  [WARN] no successful runs for {cid}")
            continue
        L = min(len(x) for x in series)
        M = np.array([x[:L] for x in series])
        per_config_weekly[cid] = M
        bw = band_width(M)
        med = np.median(M, axis=0)
        peak_h = float(med.max()); peak_w = int(med.argmax())
        attack = float(M.sum(axis=1).mean()) / (args.n_humans / 100000.0) / 1000.0
        rows.append((cid, bw, peak_h, peak_w, attack, len(series)))
        if cid == "baseline":
            base_bw, base_peak = bw, peak_h

    # Write per-config summary + NSC.
    summ = outdir / "sa_summary.csv"
    with open(summ, "w") as f:
        f.write("config,bandwidth70,peak_height,peak_week,attack_pct,n_reps,"
                "NSC_band,NSC_peak\n")
        for (cid, bw, ph, pw, atk, nr) in rows:
            nsc_b = ((bw - base_bw) / base_bw / PERTURB_DELTA) if base_bw else float("nan")
            nsc_p = ((ph - base_peak) / base_peak / PERTURB_DELTA) if base_peak else float("nan")
            f.write(f"{cid},{bw:.3f},{ph:.3f},{pw},{atk:.3f},{nr},{nsc_b:.3f},{nsc_p:.3f}\n")
    print(f"[SA] wrote {summ}")

    # Console table sorted by |NSC_band|.
    print("\n=== Sensitivity (corrected model, epidemic-regime nominal) ===")
    print(f"{'config':<12} {'band70':>8} {'peakH':>8} {'peakWk':>6} {'attack%':>8} "
          f"{'NSC_band':>9} {'NSC_peak':>9}")
    def keyf(r):
        cid, bw, ph, pw, atk, nr = r
        return abs((bw - base_bw) / base_bw) if (base_bw and cid != "baseline") else -1
    for (cid, bw, ph, pw, atk, nr) in sorted(rows, key=keyf, reverse=True):
        nsc_b = ((bw - base_bw) / base_bw / PERTURB_DELTA) if base_bw else float("nan")
        nsc_p = ((ph - base_peak) / base_peak / PERTURB_DELTA) if base_peak else float("nan")
        print(f"{cid:<12} {bw:>8.2f} {ph:>8.1f} {pw:>6} {atk:>8.2f} "
              f"{nsc_b:>9.2f} {nsc_p:>9.2f}")

    np.savez(outdir / "sa_weekly.npz", **per_config_weekly)
    print(f"[SA] wrote {outdir/'sa_weekly.npz'}")


if __name__ == "__main__":
    main()
