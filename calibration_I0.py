#!/usr/bin/env python3
"""
Calibrate I0 with batching support:
  - mode=run        : run a replicate RANGE for all (or selected) candidates
  - mode=aggregate  : aggregate existing runs, compute RMSEs, select \hat I0
  - mode=final      : run a replicate RANGE for \hat I0 only (final ensemble)

Observed CSV schema:
  week, week_continuous, incidence_rate_adjusted

One tick = one day. Weeks are 7-day non-overlapping windows starting at tick 1.
"""

from __future__ import annotations
import argparse
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd

# -----------------------------
# Observed data
# -----------------------------

def read_observed(observed_csv: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Read observed series; use week_continuous to avoid year wrap."""
    df = pd.read_csv(observed_csv)
    required = {'week', 'week_continuous', 'incidence_rate_adjusted'}
    have = {c.strip() for c in df.columns}
    missing = required.difference(have)
    if missing:
        raise ValueError(
            f"Observed CSV must contain columns {required}. "
            f"Found: {df.columns.tolist()}"
        )
    df = df.sort_values('week_continuous').reset_index(drop=True)
    weeks = df['week_continuous'].to_numpy()
    # N = 100k -> incidence per 100k equals weekly case counts
    y = df['incidence_rate_adjusted'].astype(float).to_numpy()
    return weeks, y

# -----------------------------
# Props & model invocation
# -----------------------------

def rewrite_props(props_path: Path, updates: Dict[str, str]) -> None:
    text = props_path.read_text(encoding='utf-8', errors='ignore')
    for key, val in updates.items():
        pattern = re.compile(rf"^\s*{re.escape(key)}\s*=.*$", flags=re.MULTILINE)
        if pattern.search(text):
            text = pattern.sub(f"{key}={val}", text)
        else:
            if not text.endswith("\n"):
                text += "\n"
            text += f"{key}={val}\n"
    props_path.write_text(text, encoding='utf-8')

def run_model(make_target: str = "all", workdir: Path = Path(".")) -> None:
    proc = subprocess.run(["make", make_target], cwd=str(workdir))
    if proc.returncode != 0:
        raise RuntimeError(f"`make {make_target}` failed with code {proc.returncode}")

# -----------------------------
# Output processing
# -----------------------------

def read_results_csv(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    rename_map = {}
    for col in df.columns:
        cl = col.strip().lower()
        if cl == 'tick':
            rename_map[col] = 'tick'
        if 'new' in cl and 'case' in cl and 'day' in cl:
            rename_map[col] = 'new_cases_day'
    if 'tick' not in rename_map.values():
        rename_map[df.columns[0]] = 'tick'
    if 'new_cases_day' not in rename_map.values():
        rename_map[df.columns[-1]] = 'new_cases_day'
    df = df.rename(columns=rename_map)
    if 'tick' not in df.columns or 'new_cases_day' not in df.columns:
        raise ValueError(f"results.csv missing required columns; got {df.columns.tolist()}")
    return df

def daily_to_weekly(df: pd.DataFrame, start_tick: int = 1) -> np.ndarray:
    dfd = df.sort_values('tick')
    dfd = dfd[dfd['tick'] >= start_tick].copy()
    dfd['d0'] = (dfd['tick'] - start_tick).astype(int)
    dfd['w'] = (dfd['d0'] // 7).astype(int)
    weekly = dfd.groupby('w', as_index=False)['new_cases_day'].sum().sort_values('w')
    return weekly['new_cases_day'].to_numpy()

def rmse(a: np.ndarray, b: np.ndarray) -> float:
    T = min(len(a), len(b))
    if T == 0:
        return float('inf')
    return float(np.sqrt(np.mean((a[:T] - b[:T]) ** 2)))

def peak_metrics(series: np.ndarray) -> Tuple[int, float]:
    if len(series) == 0:
        return (-1, -1.0)
    idx = int(np.argmax(series))
    val = float(series[idx])
    return idx, val

# -----------------------------
# Domain for I0
# -----------------------------

def admissible_interval(Wm1: int) -> Tuple[List[int], int, int, int]:
    LB = int(np.ceil((2.0 / 7.0) * Wm1))
    UB = int(Wm1)
    E = int(np.round((4.0 / 7.0) * Wm1))
    if UB - LB + 1 <= 5:
        domain = list(range(LB, UB + 1))
    else:
        span = 2
        lo = max(LB, E - span)
        hi = min(UB, E + span)
        while hi - lo + 1 < 5 and lo > LB:
            lo -= 1
        while hi - lo + 1 < 5 and hi < UB:
            hi += 1
        domain = list(range(lo, hi + 1))
    return domain, LB, E, UB

# -----------------------------
# Run management (batching)
# -----------------------------

def replicate_seed(base_seed: int, r: int) -> int:
    """CRN across candidates: seed depends only on r, not on I0."""
    return base_seed + r

def result_path(stash_dir: Path, I0: int, seed: int) -> Path:
    return stash_dir / f'I0_{I0:04d}' / f'seed_{seed:06d}.csv'

def run_replicate_if_needed(
    I0: int,
    r: int,
    base_seed: int,
    props_path: Path,
    make_target: str,
    workdir: Path,
    stash_dir: Path
) -> Optional[Path]:
    seed = replicate_seed(base_seed, r)
    dst = result_path(stash_dir, I0, seed)
    if dst.exists():
        # Idempotent: skip existing
        return dst
    rewrite_props(props_path, {
        'count.of.infected.humans': str(I0),
        'random.seed': str(seed),
        'global.random.seed': str(seed),  # ensure Repast RNG changes
        # 'hi.seed': str(seed)  # only if you also want HI activation to vary
    })
    run_model(make_target=make_target, workdir=workdir)
    src = workdir / 'results.csv'
    if not src.exists():
        raise FileNotFoundError(f"results.csv not found after run (I0={I0}, r={r}).")
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(src), str(dst))
    return dst

def load_candidate_weeklies(
    stash_dir: Path,
    I0: int,
    r_indices: List[int],
    base_seed: int
) -> List[np.ndarray]:
    weeklies = []
    for r in r_indices:
        seed = replicate_seed(base_seed, r)
        p = result_path(stash_dir, I0, seed)
        if not p.exists():
            continue
        df = read_results_csv(p)
        weeklies.append(daily_to_weekly(df, start_tick=1))
    return weeklies

# -----------------------------
# Modes
# -----------------------------

def mode_run(args, domain: List[int]) -> None:
    """Run a batch: replicates in [r_start, r_end) for given domain."""
    r_indices = list(range(args.r_start, args.r_end))
    print(f"[RUN] Candidates={domain}, replicates={r_indices}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    for I0 in domain:
        for r in r_indices:
            try:
                path = run_replicate_if_needed(
                    I0=I0,
                    r=r,
                    base_seed=args.base_seed,
                    props_path=args.props,
                    make_target=args.make_target,
                    workdir=args.workdir,
                    stash_dir=args.outdir,
                )
                if path is not None:
                    print(f"  I0={I0} r={r} -> {path}")
                else:
                    print(f"  I0={I0} r={r} -> exists (skipped)")
            except Exception as e:
                print(f"  I0={I0} r={r} -> ERROR: {e}")

def mode_aggregate(args, domain: List[int], weeks_eval: np.ndarray, y_eval: np.ndarray) -> None:
    """Aggregate whatever is on disk, compute RMSEs, and select \hat I0."""
    print("[AGGREGATE] Reading existing runs from disk...")
    records: List[Dict[str, object]] = []
    # use full intended replicate set R for comparability (ignore missing ones)
    intended_r = list(range(args.replicates))
    for I0 in domain:
        weeklies = load_candidate_weeklies(args.outdir, I0, intended_r, args.base_seed)
        if not weeklies:
            print(f"  I0={I0}: no data found; skipping.")
            continue
        T_star = min(len(w) for w in weeklies)
        W = np.vstack([w[:T_star] for w in weeklies])
        mean_w = W.mean(axis=0)
        T = min(len(mean_w), len(y_eval))
        J = rmse(mean_w[:T], y_eval[:T])
        p_idx, p_val = peak_metrics(mean_w[:T])
        records.append({
            'I0': I0, 'J': J, 'peak_week_index': p_idx, 'peak_value': p_val,
            'T_used': T, 'weekly_mean': mean_w[:T]
        })
        print(f"  I0={I0}: replicates={len(weeklies)} RMSE={J:.3f}")

    if not records:
        raise RuntimeError("No candidate data to aggregate.")

    # Tie-breaks
    def tie_key(rec):
        T_used = rec['T_used']
        if T_used > 0:
            peak_idx_obs = int(np.argmax(y_eval[:T_used]))
            peak_val_obs = float(np.max(y_eval[:T_used]))
        else:
            peak_idx_obs = -1
            peak_val_obs = -1.0
        return (
            float(rec['J']),
            abs(int(rec['peak_week_index']) - peak_idx_obs),
            abs(float(rec['peak_value']) - peak_val_obs)
        )

    records.sort(key=tie_key)
    best = records[0]
    hat_I0 = int(best['I0'])
    print(f"\n==> Selected I0 = {hat_I0}  RMSE = {best['J']:.3f}")

    # Write summary table
    with open(args.outdir / 'selection_summary.csv', 'w', encoding='utf-8') as fh:
        fh.write('I0,J,peak_week_index,peak_value,T_used\n')
        for rec in records:
            fh.write(f"{rec['I0']},{float(rec['J']):.6f},{int(rec['peak_week_index'])},"
                     f"{float(rec['peak_value']):.6f},{int(rec['T_used'])}\n")

    # Persist choice for final mode convenience
    (args.outdir / 'hat_I0.txt').write_text(str(hat_I0), encoding='utf-8')

def mode_final(args, weeks_eval: np.ndarray, y_eval: np.ndarray) -> None:
    """Run final ensemble for \hat I0 in batches, then (optionally) aggregate when complete."""
    # Determine \hat I0
    if args.hat_I0 is not None:
        hat = int(args.hat_I0)
    else:
        hat_file = args.outdir / 'hat_I0.txt'
        if not hat_file.exists():
            raise RuntimeError("hat_I0 not provided and hat_I0.txt not found. Run --mode aggregate first or pass --hat-I0.")
        hat = int(hat_file.read_text().strip())
    print(f"[FINAL] I0={hat}  replicates [{args.r_start}, {args.r_end}) (of R_final={args.replicates_final})")

    # Run the requested batch for final ensemble (distinct seed block)
    r_indices = list(range(args.r_start, args.r_end))
    final_base = args.base_seed + 10_000
    args.outdir.mkdir(parents=True, exist_ok=True)

    for r in r_indices:
        try:
            path = run_replicate_if_needed(
                I0=hat, r=r, base_seed=final_base,
                props_path=args.props, make_target=args.make_target,
                workdir=args.workdir, stash_dir=args.outdir
            )
            if path is not None:
                print(f"  I0={hat} r={r} -> {path}")
            else:
                print(f"  I0={hat} r={r} -> exists (skipped)")
        except Exception as e:
            print(f"  I0={hat} r={r} -> ERROR: {e}")

    # If user requests aggregation now (only if full set present)
    have_all = all(result_path(args.outdir, hat, replicate_seed(final_base, r)).exists()
                   for r in range(args.replicates_final))
    if have_all:
        print("[FINAL] Aggregating full final ensemble...")
        weeklies = load_candidate_weeklies(args.outdir, hat, list(range(args.replicates_final)), final_base)
        T_star = min(len(w) for w in weeklies)
        W = np.vstack([w[:T_star] for w in weeklies])
        T = min(T_star, len(y_eval))
        out = pd.DataFrame({
            'week_idx': np.arange(T),
            'obs': y_eval[:T],
            'median': np.median(W, axis=0)[:T],
            'p15': np.percentile(W, 15, axis=0)[:T],
            'p85': np.percentile(W, 85, axis=0)[:T],
            'mean': W.mean(axis=0)[:T],
        })
        out_csv = args.outdir / f'final_ensemble_I0_{hat:04d}.csv'
        out.to_csv(out_csv, index=False)
        print(f"[FINAL] Wrote {out_csv}")
    else:
        print("[FINAL] Not aggregating yet (missing some replicates).")

# -----------------------------
# Main
# -----------------------------

def main():
    ap = argparse.ArgumentParser(description="Calibrate I0 with batching.")
    ap.add_argument('--mode', choices=['run','aggregate','final'], required=True)
    ap.add_argument('--observed', required=True, type=Path,
                    help='Path to observed CSV (week, week_continuous, incidence_rate_adjusted).')
    ap.add_argument('--props', required=True, type=Path,
                    help='Path to props/model.props to rewrite.')
    ap.add_argument('--config', required=True, type=Path,
                    help='Path to props/config.props (not modified, but required by your workflow).')
    ap.add_argument('--workdir', type=Path, default=Path('.'),
                    help='Project root where `make` is executed.')
    ap.add_argument('--make-target', type=str, default='all',
                    help='Make target that runs the simulation.')
    ap.add_argument('--replicates', type=int, default=32,
                    help='R: replicates per candidate during SELECTION.')
    ap.add_argument('--replicates-final', type=int, default=64,
                    help='R_final: replicates for FINAL ensemble.')
    ap.add_argument('--base-seed', type=int, default=12345,
                    help='Base seed for CRN.')
    ap.add_argument('--outdir', type=Path, default=Path('calib_I0_runs'),
                    help='Directory to store run outputs.')
    ap.add_argument('--r-start', type=int, default=0,
                    help='Replicate start index (inclusive) for this batch.')
    ap.add_argument('--r-end', type=int, default=0,
                    help='Replicate end index (exclusive) for this batch.')
    ap.add_argument('--week-min', type=int, default=None,
                    help='Optional: minimum continuous week to include in RMSE window.')
    ap.add_argument('--week-max', type=int, default=None,
                    help='Optional: maximum continuous week to include in RMSE window.')
    ap.add_argument('--candidates', type=str, default=None,
                    help='Optional comma-separated I0 list to run in mode=run (e.g., "30,31,32").')
    ap.add_argument('--hat-I0', type=int, default=None,
                    help='Optional explicit \u02C6I0 for mode=final; otherwise read from hat_I0.txt.')
    args = ap.parse_args()

    # Read observed series and fix evaluation window
    weeks_obs, y_all = read_observed(args.observed)
    if args.week_min is not None or args.week_max is not None:
        wmin = args.week_min if args.week_min is not None else int(np.min(weeks_obs))
        wmax = args.week_max if args.week_max is not None else int(np.max(weeks_obs))
        mask = (weeks_obs >= wmin) & (weeks_obs <= wmax)
        weeks = weeks_obs[mask]
        y = y_all[mask]
    else:
        weeks, y = weeks_obs, y_all

    # Determine W_-1 (continuous week 26 if present; else first point)
    if np.min(weeks_obs) <= 26 and np.any(weeks_obs == 26):
        Wm1 = int(round(float(y_all[np.where(weeks_obs == 26)[0][0]])))
        mask_eval = weeks_obs > 26
        y_eval = y_all[mask_eval]
        weeks_eval = weeks_obs[mask_eval]
    else:
        Wm1 = int(round(float(y[0])))
        y_eval = y
        weeks_eval = weeks

    # Domain for selection
    domain, LB, E, UB = admissible_interval(Wm1)
    if args.candidates is not None:
        # Allow user to restrict domain manually
        dom_user = [int(x.strip()) for x in args.candidates.split(',') if x.strip()]
        domain = dom_user

    print(f"Domain={domain}  LB={LB} E={E} UB={UB}")

    if args.mode == 'run':
        if args.r_end <= args.r_start:
            raise ValueError("Provide a positive replicate range with --r-start and --r-end.")
        mode_run(args, domain)

    elif args.mode == 'aggregate':
        mode_aggregate(args, domain, weeks_eval, y_eval)

    elif args.mode == 'final':
        if args.r_end <= args.r_start:
            raise ValueError("Provide a positive replicate range with --r-start and --r-end.")
        mode_final(args, weeks_eval, y_eval)

if __name__ == '__main__':
    main()
