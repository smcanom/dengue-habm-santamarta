# Quick Start Guide: Sensitivity Analysis

## Prerequisites

- ✅ Repast HPC installed
- ✅ Model compiles successfully
- ✅ Python 3.6+ with pandas and numpy

## 5-Minute Test

### 1. Verify Parameter Parsing

```bash
python3 test_params.py
```

**Expected output**:
```
✓ All parameter tests passed!
```

### 2. Compile Model

```bash
make clean
make
```

**Expected**: No compilation errors

### 3. Run Minimal Sensitivity Test (5 runs)

```bash
# 5 configs × 1 replicate = 5 runs (~5-10 minutes)
python3 run_sensitivity.py --mode run --replicates 1 --outdir test_quick
```

**Expected**: 5 successful runs in `test_quick/`

### 4. Aggregate Results

```bash
python3 run_sensitivity.py --mode aggregate --outdir test_quick --replicates 1
```

**Expected output files**:
- `test_quick/sensitivity_indices.csv`
- `test_quick/config_pi70.csv`
- `test_quick/all_weekly_cases.csv`

### 5. View Results

```bash
cat test_quick/sensitivity_indices.csv
```

**Expected columns**: `param_name`, `NSC`, `W_baseline`, `W_perturbed`, `delta_rmse`, `delta_cov`

---

## Full Analysis (Recommended)

### Run with 3 Replicates

```bash
# 5 configs × 3 replicates = 15 runs (~30-60 minutes)
python3 run_sensitivity.py --mode run --replicates 3 --outdir out_sensitivity

# Aggregate
python3 run_sensitivity.py --mode aggregate --outdir out_sensitivity --replicates 3

# View
cat out_sensitivity/sensitivity_indices.csv
```

### Interpret Results

**Sensitivity Index (NSC)**:
- `NSC > 0.5`: **High influence** on prediction uncertainty
- `0.1 < NSC < 0.5`: **Moderate influence**
- `NSC < 0.1`: **Low influence**
- `NSC < 0`: Parameter **reduces** uncertainty (stabilizes model)

**Example output**:
```
param_name,NSC,W_baseline,W_perturbed,delta_rmse,delta_cov
sigma_M,0.25,5.234,5.891,0.757,0.0
sigma_H,0.41,5.234,6.389,1.234,0.0
z,0.19,5.234,5.821,0.512,0.0
beta_mh,0.57,5.234,6.719,1.892,0.0
```

**Interpretation**:
- `beta_mh` has the **highest sensitivity** (NSC=0.57)
- 10% increase in `beta_mh` → 57% increase in uncertainty width
- `sigma_M` has the **lowest sensitivity** (NSC=0.25)

---

## Troubleshooting

### "Command 'make' not found"

**Solution**: Install build tools
```bash
sudo apt-get install build-essential
```

### "Module 'pandas' not found"

**Solution**: Install Python dependencies
```bash
pip3 install pandas numpy --user
```

### "Model exited with code 1"

**Solution**: Check simulation parameters
```bash
# Reduce population for faster testing
nano props/model.props
# Change: count.of.humans = 10000
```

### "weekly_neighborhood_cases.csv not found"

**Solution**: Increase simulation time
```bash
nano props/model.props
# Change: stop.at = 365
```

---

## Command Reference

### Test Mode (1 replicate)
```bash
python3 run_sensitivity.py --mode run --replicates 1 --outdir test_run
python3 run_sensitivity.py --mode aggregate --outdir test_run --replicates 1
```

### Quick Mode (3 replicates)
```bash
python3 run_sensitivity.py --mode run --replicates 3
python3 run_sensitivity.py --mode aggregate --replicates 3
```

### Full Mode (20 replicates)
```bash
python3 run_sensitivity.py --mode run --replicates 20 --outdir out_full
python3 run_sensitivity.py --mode aggregate --outdir out_full --replicates 20
```

### Custom Parameters
```bash
python3 run_sensitivity.py --mode run \
    --replicates 5 \
    --model_bin ./Main.exe \
    --config_props props/config.props \
    --model_props props/model.props \
    --outdir my_analysis \
    --obs_csv view/incidence_rates_santa_marta.csv
```

---

## File Structure After Running

```
out_sensitivity/
├── baseline/
│   ├── replicate_0/
│   │   ├── model.props
│   │   ├── out/
│   │   │   ├── weekly_cases.csv
│   │   │   └── run_summary.csv
│   │   ├── weekly_neighborhood_cases.csv
│   │   └── results.csv
│   ├── replicate_1/
│   └── replicate_2/
├── sigma_M_p10/
│   ├── replicate_0/
│   ├── replicate_1/
│   └── replicate_2/
├── sigma_H_p10/
├── z_p10/
├── beta_mh_p10/
├── all_weekly_cases.csv          ← Combined time series
├── all_run_summaries.csv         ← RMSE summaries
├── config_pi70.csv                ← Prediction intervals
├── config_fit.csv                 ← Aggregate metrics
├── sensitivity_indices.csv        ← Main results
└── meta.json                      ← Metadata
```

---

## Next Steps

1. ✅ **Quick test**: `python3 test_params.py`
2. ✅ **Compile**: `make clean && make`
3. ✅ **Test run**: `python3 run_sensitivity.py --mode run --replicates 1 --outdir test_run`
4. ✅ **Full analysis**: `python3 run_sensitivity.py --mode run --replicates 3`
5. 📊 **Visualize**: See `SENSITIVITY_ANALYSIS_README.md` for plotting examples

---

## Help

- **Full documentation**: `SENSITIVITY_ANALYSIS_README.md`
- **Implementation details**: `IMPLEMENTATION_SUMMARY.md`
- **Original workflow reference**: `calibration_I0.py`

**Need help?** Check the READMEs or inspect `out_sensitivity/meta.json` for run details.

