#!/usr/bin/env python3
"""
run_sensitivity.py – Minimal SA/UA pipeline for Repast HPC dengue model

One-sided OAT perturbations with CRN:
  - Baseline + four +10% perturbations
  - R replicates per configuration
  - Common Random Numbers (replicate_id = base_seed + r)

Usage:
  python run_sensitivity.py --mode run --r-start 0 --r-end 1
  python run_sensitivity.py --mode aggregate

Pattern follows calibration_I0.py:
  1. Edit props/model.props in place
  2. Run 'make all' from project root
  3. Move results.csv to replicate folder
"""

from __future__ import annotations
import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd


# ============================================================================
# Configuration
# ============================================================================

NOMINAL_PARAMS = {
    'sigma_M': 0.3,
    'sigma_H': 3.0,
    'z': 0.3,
    'r': 0.6,
    'C': 30.0,
    'beta_mh': 0.10,
    'beta_hm': 0.10,
}

PERTURB_PARAMS = ['sigma_M', 'sigma_H', 'z', 'beta_mh']  # OAT
PERTURB_DELTA = 0.10  # +10%

BASE_SEED = 12345

# Observed data for RMSE
OBS_CSV = 'view/incidence_rates_santa_marta.csv'  # relative to project root


# ============================================================================
# Utilities (following calibration_I0.py pattern)
# ============================================================================

def rewrite_props(props_path: Path, updates: Dict[str, str]) -> None:
    """Edit key=value pairs in a properties file in place."""
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
    """Run 'make <target>' from workdir. Raises on failure."""
    proc = subprocess.run(["make", make_target], cwd=str(workdir))
    if proc.returncode != 0:
        raise RuntimeError(f"`make {make_target}` failed with code {proc.returncode}")


def replicate_seed(base_seed: int, r: int) -> int:
    """CRN: seed depends only on replicate_id, not on configuration."""
    return base_seed + r


def result_path(stash_dir: Path, config_id: str, seed: int) -> Path:
    """Path where results.csv for (config_id, seed) will be stored."""
    return stash_dir / f'{config_id}' / f'seed_{seed:06d}.csv'


# ============================================================================
# Configuration Generation
# ============================================================================

def generate_configs() -> List[Dict[str, object]]:
    """
    Build configuration list:
      - baseline
      - one +10% perturbation per PERTURB_PARAMS
    """
    configs = []

    # Baseline
    configs.append({
        'config_id': 'baseline',
        'perturb_param': 'baseline',
        'perturb_delta': 0.0,
        'params': NOMINAL_PARAMS.copy(),
    })

    # OAT perturbations (one-sided)
    for p in PERTURB_PARAMS:
        perturbed = NOMINAL_PARAMS.copy()
        perturbed[p] = NOMINAL_PARAMS[p] * (1.0 + PERTURB_DELTA)
        configs.append({
            'config_id': f'{p}_p10',
            'perturb_param': p,
            'perturb_delta': PERTURB_DELTA,
            'params': perturbed,
        })

    return configs


# ============================================================================
# Run Management
# ============================================================================

def run_replicate_if_needed(
    config_id: str,
    params: Dict[str, float],
    perturb_param: str,
    perturb_delta: float,
    r: int,
    base_seed: int,
    props_path: Path,
    make_target: str,
    workdir: Path,
    stash_dir: Path,
    obs_csv: str
) -> Optional[Path]:
    """
    Run one replicate for a configuration.
    1. Edit props/model.props in place
    2. Run 'make all' from project root
    3. Move results.csv → stash_dir/<config_id>/seed_<seed>.csv
    
    Returns destination path if run succeeded, None if skipped (already exists).
    """
    seed = replicate_seed(base_seed, r)
    dst = result_path(stash_dir, config_id, seed)
    
    if dst.exists():
        # Idempotent: skip existing
        return None
    
    # Prepare model.props updates
    updates = {}
    for k, v in params.items():
        updates[k] = str(v)
    
    updates['base_seed'] = str(base_seed)
    updates['replicate_id'] = str(r)
    updates['config_id'] = config_id
    updates['perturb_param'] = perturb_param
    updates['perturb_delta'] = str(perturb_delta)
    updates['obs_csv'] = obs_csv
    
    # Also set random.seed for Repast RNG (CRN)
    updates['random.seed'] = str(seed)
    updates['global.random.seed'] = str(seed)
    
    # Edit props/model.props in place
    rewrite_props(props_path, updates)
    
    # Run model from project root (all input data accessed naturally)
    run_model(make_target=make_target, workdir=workdir)
    
    # Move results.csv to stash
    src = workdir / 'results.csv'
    if not src.exists():
        raise FileNotFoundError(f"results.csv not found after run (config={config_id}, r={r}).")
    
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(src), str(dst))
    
    return dst


# ============================================================================
# Mode: run
# ============================================================================

def mode_run(args) -> None:
    """
    Run a batch: replicates in [r_start, r_end) for all configurations.
    """
    configs = generate_configs()
    r_indices = list(range(args.r_start, args.r_end))
    
    print(f"[RUN] Configs={len(configs)}, replicates={r_indices}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    
    for cfg in configs:
        for r in r_indices:
            try:
                path = run_replicate_if_needed(
                    config_id=cfg['config_id'],
                    params=cfg['params'],
                    perturb_param=cfg['perturb_param'],
                    perturb_delta=cfg['perturb_delta'],
                    r=r,
                    base_seed=args.base_seed,
                    props_path=args.props,
                    make_target=args.make_target,
                    workdir=args.workdir,
                    stash_dir=args.outdir,
                    obs_csv=args.obs_csv,
                )
                if path is not None:
                    print(f"  {cfg['config_id']} r={r} -> {path}")
                else:
                    print(f"  {cfg['config_id']} r={r} -> exists (skipped)")
            except Exception as e:
                print(f"  {cfg['config_id']} r={r} -> ERROR: {e}")


# ============================================================================
# Output Processing (for aggregation)
# ============================================================================

def read_results_csv(csv_path: Path) -> pd.DataFrame:
    """Read and normalize results.csv columns."""
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
    """Convert daily case counts to weekly totals."""
    dfd = df.sort_values('tick')
    dfd = dfd[dfd['tick'] >= start_tick].copy()
    dfd['d0'] = (dfd['tick'] - start_tick).astype(int)
    dfd['w'] = (dfd['d0'] // 7).astype(int)
    weekly = dfd.groupby('w', as_index=False)['new_cases_day'].sum().sort_values('w')
    return weekly['new_cases_day'].to_numpy()


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    """Root mean squared error between two series."""
    T = min(len(a), len(b))
    if T == 0:
        return float('inf')
    return float(np.sqrt(np.mean((a[:T] - b[:T]) ** 2)))


def load_config_weeklies(
    stash_dir: Path,
    config_id: str,
    r_indices: List[int],
    base_seed: int
) -> List[np.ndarray]:
    """Load weekly time series for all replicates of a configuration."""
    weeklies = []
    for r in r_indices:
        seed = replicate_seed(base_seed, r)
        p = result_path(stash_dir, config_id, seed)
        if not p.exists():
            continue
        df = read_results_csv(p)
        weeklies.append(daily_to_weekly(df, start_tick=1))
    return weeklies


# ============================================================================
# Mode: aggregate
# ============================================================================

def mode_aggregate(args) -> None:
    """
    Aggregate all existing runs:
      - Compute per-config PI bands, RMSE, coverage, W
      - Compute sensitivity indices (NSC, delta_rmse, delta_cov)
    """
    configs = generate_configs()
    print("[AGGREGATE] Aggregating results...")
    
    # Load observed data
    obs_df = None
    y_obs = None
    weeks_obs = None
    if Path(args.obs_csv).exists():
        obs_df = pd.read_csv(args.obs_csv)
        weeks_obs = obs_df['week'].values if 'week' in obs_df.columns else obs_df.iloc[:, 0].values
        y_obs = obs_df['incidence_rate_adjusted'].values if 'incidence_rate_adjusted' in obs_df.columns else obs_df.iloc[:, 1].values
    
    # Intended replicate set (may not all be present)
    intended_r = list(range(args.replicates))
    
    config_fit = []
    all_weekly_records = []
    
    for cfg in configs:
        cid = cfg['config_id']
        weeklies = load_config_weeklies(args.outdir, cid, intended_r, args.base_seed)
        
        if not weeklies:
            print(f"  {cid}: no data found; skipping.")
            continue
        
        # Align all replicates to shortest length
        T_star = min(len(w) for w in weeklies)
        W = np.vstack([w[:T_star] for w in weeklies])
        
        mean_w = W.mean(axis=0)
        p15 = np.percentile(W, 15, axis=0)
        p85 = np.percentile(W, 85, axis=0)
        
        # Mean prediction interval width
        mean_W = float(np.mean(p85 - p15))
        
        # RMSE against observed (if available)
        mean_rmse = np.nan
        if y_obs is not None:
            T = min(len(mean_w), len(y_obs))
            mean_rmse = rmse(mean_w[:T], y_obs[:T])
        
        # Coverage70 (fraction of observed points inside [p15, p85])
        coverage70 = np.nan
        if y_obs is not None and weeks_obs is not None:
            common_weeks = np.arange(min(len(mean_w), len(y_obs)))
            if len(common_weeks) > 0:
                in_band = sum(1 for i in common_weeks if p15[i] <= y_obs[i] <= p85[i])
                coverage70 = in_band / len(common_weeks)
        
        config_fit.append({
            'config_id': cid,
            'mean_rmse': mean_rmse,
            'Coverage70': coverage70,
            'mean_W': mean_W,
            'n_replicates': len(weeklies),
        })
        
        # Store weekly time series for all_weekly_cases.csv
        for r_idx, r in enumerate(intended_r):
            if r_idx < len(weeklies):
                for week_idx, cases in enumerate(weeklies[r_idx]):
                    all_weekly_records.append({
                        'config_id': cid,
                        'replicate_id': r,
                        'week': week_idx,
                        'cases': cases,
                    })
        
        print(f"  {cid}: replicates={len(weeklies)}, RMSE={mean_rmse:.3f}, W={mean_W:.2f}, Cov70={coverage70:.3f}")
    
    if not config_fit:
        raise RuntimeError("No configuration data to aggregate.")
    
    # Write per-config fit summary
    df_config_fit = pd.DataFrame(config_fit)
    fit_path = args.outdir / 'config_fit.csv'
    df_config_fit.to_csv(fit_path, index=False)
    print(f"  Wrote {fit_path}")
    
    # Write all weekly cases
    df_all_weekly = pd.DataFrame(all_weekly_records)
    weekly_path = args.outdir / 'all_weekly_cases.csv'
    df_all_weekly.to_csv(weekly_path, index=False)
    print(f"  Wrote {weekly_path}")
    
    # Sensitivity indices
    baseline_row = df_config_fit[df_config_fit['config_id'] == 'baseline']
    if not baseline_row.empty:
        W_b = baseline_row['mean_W'].values[0]
        rmse_b = baseline_row['mean_rmse'].values[0]
        cov_b = baseline_row['Coverage70'].values[0]
        
        sens_records = []
        for param in PERTURB_PARAMS:
            perturb_id = f'{param}_p10'
            perturb_row = df_config_fit[df_config_fit['config_id'] == perturb_id]
            if not perturb_row.empty:
                W_p = perturb_row['mean_W'].values[0]
                rmse_p = perturb_row['mean_rmse'].values[0]
                cov_p = perturb_row['Coverage70'].values[0]
                
                # NSC = (W_p - W_b) / W_b / delta
                nsc = ((W_p - W_b) / W_b) / PERTURB_DELTA if W_b > 0 else np.nan
                delta_rmse = rmse_p - rmse_b
                delta_cov = cov_p - cov_b if not np.isnan(cov_p) and not np.isnan(cov_b) else np.nan
                
                sens_records.append({
                    'parameter': param,
                    'NSC': nsc,
                    'W_baseline': W_b,
                    'W_perturbed': W_p,
                    'delta_rmse': delta_rmse,
                    'delta_cov': delta_cov,
                })
        
        if sens_records:
            df_sens = pd.DataFrame(sens_records)
            sens_path = args.outdir / 'sensitivity_indices.csv'
            df_sens.to_csv(sens_path, index=False)
            print(f"  Wrote {sens_path}")
    
    # Write meta.json
    meta = {
        'base_seed': args.base_seed,
        'replicates': args.replicates,
        'nominal_params': NOMINAL_PARAMS,
        'perturb_params': PERTURB_PARAMS,
        'perturb_delta': PERTURB_DELTA,
        'obs_csv': args.obs_csv,
    }
    meta_path = args.outdir / 'meta.json'
    with open(meta_path, 'w') as f:
        json.dump(meta, f, indent=2)
    print(f"  Wrote {meta_path}")
    
    print("[AGGREGATE] Done.")


# ============================================================================
# Main
# ============================================================================

def main():
    ap = argparse.ArgumentParser(description="Minimal SA/UA pipeline for dengue model.")
    ap.add_argument('--mode', choices=['run', 'aggregate'], required=True)
    ap.add_argument('--props', type=Path, default=Path('props/model.props'),
                    help='Path to model.props to rewrite in place.')
    ap.add_argument('--workdir', type=Path, default=Path('.'),
                    help='Project root where `make` is executed.')
    ap.add_argument('--make-target', type=str, default='all',
                    help='Make target that runs the simulation.')
    ap.add_argument('--base-seed', type=int, default=BASE_SEED,
                    help='Base seed for CRN.')
    ap.add_argument('--replicates', type=int, default=1,
                    help='Total number of replicates per configuration.')
    ap.add_argument('--outdir', type=Path, default=Path('test_run'),
                    help='Directory to store run outputs.')
    ap.add_argument('--r-start', type=int, default=0,
                    help='Replicate start index (inclusive) for this batch.')
    ap.add_argument('--r-end', type=int, default=1,
                    help='Replicate end index (exclusive) for this batch.')
    ap.add_argument('--obs-csv', type=str, default=OBS_CSV,
                    help='Path to observed CSV for RMSE computation.')
    args = ap.parse_args()
    
    if args.mode == 'run':
        if args.r_end <= args.r_start:
            raise ValueError("Provide a positive replicate range with --r-start and --r-end.")
        mode_run(args)
    
    elif args.mode == 'aggregate':
        mode_aggregate(args)


if __name__ == '__main__':
    main()
