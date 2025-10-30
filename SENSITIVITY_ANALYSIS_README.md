# Sensitivity Analysis & Uncertainty Quantification Pipeline

## Overview

This implementation enables **parameter-driven sensitivity analysis** with **common random numbers (CRN)** for the dengue HABM model. The pipeline supports:

- **One-sided One-At-a-Time (OAT)** sensitivity analysis with +10% perturbations
- **Common Random Numbers (CRN)** for variance reduction across parameter configurations
- **Budget-aware design** focusing on 4 key parameters: `sigma_M`, `sigma_H`, `z`, `beta_mh`
- **Analysis-ready CSV outputs** for visualization and further analysis

## Architecture

### C++ Changes

#### 1. Parameter Carrier (`Params` struct)

A lightweight struct in `MyModel.h` holds all sensitivity-relevant parameters:

```cpp
struct Params {
    double sigma_M{0.3};      // mosquito mortality sensitivity
    double sigma_H{3.0};      // human bite sensitivity
    double z{0.3};            // bite aggregation parameter
    double r{0.6};            // reproduction rate
    double C{30.0};           // carrying capacity
    double beta_mh{0.10};     // mosquito-to-human transmission
    double beta_hm{0.10};     // human-to-mosquito transmission
    int base_seed{12345};     // base random seed
    int replicate_id{0};      // replicate identifier (for CRN)
    std::string config_id{"baseline"};
    std::string perturb_param{"baseline"};
    double perturb_delta{0.0};
    std::string obs_csv{""};  // path to observed data
};
```

#### 2. Parameter Threading

- **Main.cpp**: Reads parameters from `model.props` and applies perturbations before constructing the model
- **MyModel**: Stores `params_` and passes it to `SEIModel` and `Human` instances
- **SEIModel**: Uses `params_->r`, `params_->C`, `params_->z`, `params_->beta_hm` instead of hardcoded constants
- **Human**: Uses `params_->beta_mh` for infection probability calculations

#### 3. Random Seed Management (CRN)

The run seed is computed as:
```cpp
run_seed = params_.base_seed + params_.replicate_id
```

**Critical**: The same `replicate_id` values (e.g., 0, 1, 2) are used across ALL configurations to enable CRN. This ensures that random number streams are synchronized across baseline and perturbed runs, reducing variance in sensitivity estimates.

#### 4. Output Files

At the end of each simulation, the model writes:

- **`out/weekly_cases.csv`**: City-wide weekly case counts with metadata
  ```
  config_id,param_name,delta_sign,replicate_id,week,cases
  baseline,baseline,0,0,1,15.2
  baseline,baseline,0,0,2,18.7
  ...
  sigma_M_p10,sigma_M,1,0,1,16.3
  ...
  ```

- **`out/run_summary.csv`**: Per-replicate summary statistics
  ```
  config_id,replicate_id,rmse,n_weeks,seed
  baseline,0,12.345678,52,12345
  baseline,1,11.987654,52,12346
  ...
  ```

## Python Driver Script

### Usage

#### Step 1: Run Simulations

```bash
python run_sensitivity.py --mode run \
    --replicates 3 \
    --model_bin ./Main.exe \
    --config_props props/config.props \
    --model_props props/model.props \
    --outdir out_sensitivity \
    --obs_csv view/incidence_rates_santa_marta.csv
```

This will:
1. Create 5 configurations (baseline + 4 perturbations)
2. Run 3 replicates for each configuration (15 runs total)
3. Store outputs in `out_sensitivity/{config_id}/replicate_{id}/`

**Configuration Structure:**
```
out_sensitivity/
├── baseline/
│   ├── replicate_0/
│   │   ├── model.props
│   │   ├── out/
│   │   │   ├── weekly_cases.csv
│   │   │   └── run_summary.csv
│   │   └── ...
│   ├── replicate_1/
│   └── replicate_2/
├── sigma_M_p10/
│   ├── replicate_0/
│   ├── replicate_1/
│   └── replicate_2/
├── sigma_H_p10/
├── z_p10/
└── beta_mh_p10/
```

#### Step 2: Aggregate Results

```bash
python run_sensitivity.py --mode aggregate \
    --outdir out_sensitivity \
    --replicates 3
```

This produces analysis-ready CSVs:

1. **`all_weekly_cases.csv`**: Combined weekly data from all runs
2. **`all_run_summaries.csv`**: Combined summary statistics
3. **`config_pi70.csv`**: Per-configuration 70% prediction intervals
   ```
   config_id,week,q15,median,q85,mean,pi_width
   baseline,1,12.5,15.2,18.1,15.3,5.6
   ...
   ```
4. **`config_fit.csv`**: Per-configuration aggregate metrics
   ```
   config_id,W,RMSE_mean,Coverage70
   baseline,5.234,12.345,-1.0
   sigma_M_p10,5.891,13.102,-1.0
   ...
   ```
5. **`sensitivity_indices.csv`**: Normalized sensitivity coefficients
   ```
   param_name,NSC,W_baseline,W_perturbed,delta_rmse,delta_cov
   sigma_M,0.251,5.234,5.891,0.757,0.0
   sigma_H,0.412,5.234,6.389,1.234,0.0
   z,0.189,5.234,5.821,0.512,0.0
   beta_mh,0.567,5.234,6.719,1.892,0.0
   ```
6. **`meta.json`**: Metadata for reproducibility

### Sensitivity Index (NSC)

The **Normalized Sensitivity Coefficient** measures the relative change in uncertainty (PI width) per unit change in parameter:

```
NSC_i = (W(θ_i^+) - W(θ_0)) / (δ × W(θ_0))
```

Where:
- `W(θ)` = mean 70% prediction interval width across weeks
- `θ_i^+` = parameter i increased by δ = 0.10 (10%)
- `θ_0` = baseline parameter set

**Interpretation:**
- `NSC > 0`: Increasing the parameter increases uncertainty
- `NSC < 0`: Increasing the parameter decreases uncertainty (stabilizes predictions)
- `|NSC|` large: Parameter is highly influential on prediction uncertainty

## Model Parameters Reference

| Parameter | Symbol | Description | Default | Used In |
|-----------|--------|-------------|---------|---------|
| `sigma_M` | σ_M | Mosquito mortality sensitivity | 0.3 | (Placeholder for future use) |
| `sigma_H` | σ_H | Human bite sensitivity | 3.0 | (Placeholder for future use) |
| `z` | z | Bite aggregation parameter | 0.3 | SEIModel (human-mosquito interaction) |
| `r` | r | Mosquito reproduction rate | 0.6 | SEIModel (egg laying) |
| `C` | C | Mosquito carrying capacity | 30.0 | SEIModel (density dependence) |
| `beta_mh` | β_mh | Mosquito→Human transmission prob. | 0.10 | Human (infection calculation) |
| `beta_hm` | β_hm | Human→Mosquito transmission prob. | 0.10 | SEIModel (new infections) |

**Note**: Parameters `sigma_M` and `sigma_H` are included in the infrastructure but not yet mapped to specific model equations. They can be connected to mortality or bite rate modifiers in future iterations.

## Workflow Example

### Quick Start (3 replicates)

```bash
# 1. Compile the model
make clean
make

# 2. Run sensitivity analysis (5 configs × 3 reps = 15 runs)
python run_sensitivity.py --mode run --replicates 3

# 3. Aggregate results
python run_sensitivity.py --mode aggregate

# 4. Explore outputs
ls out_sensitivity/
cat out_sensitivity/sensitivity_indices.csv
```

### Full Analysis (20 replicates)

For publication-quality results with tighter confidence intervals:

```bash
# Run with more replicates (may take several hours)
python run_sensitivity.py --mode run --replicates 20 \
    --outdir out_sensitivity_full

# Aggregate
python run_sensitivity.py --mode aggregate \
    --outdir out_sensitivity_full \
    --replicates 20
```

### Custom Perturbations

To test different perturbation levels, edit `run_sensitivity.py`:

```python
# Line 35
DELTA = 0.20  # Test +20% instead of +10%
```

Or add two-sided perturbations:

```python
# In mode_run(), add negative deltas:
configs.append({
    "config_id": f"{param}_m10",
    "perturb_param": param,
    "perturb_delta": -0.10  # -10%
})
```

## Integration with Existing Workflow

### Compatibility with `calibration_I0.py`

The sensitivity pipeline is **independent** of `calibration_I0.py` but follows similar design patterns:

- Both use `model.props` for parameter control
- Both support batching via command-line arguments
- Both write analysis-ready CSVs

You can:
1. Run `calibration_I0.py` to find optimal `I0`
2. Update `count.of.infected.humans` in `model.props`
3. Run `run_sensitivity.py` with the calibrated initial conditions

### Custom Observed Data

If you have alternative observed data (e.g., different time windows or neighborhoods):

1. Create a CSV with columns: `week`, `week_continuous`, `incidence_rate_adjusted`
2. Pass it via `--obs_csv path/to/your_data.csv`
3. The C++ code will compute RMSE against this file

## Output Visualization (Future Work)

The CSVs are ready for plotting. Suggested visualizations:

1. **Tornado Diagram**: Plot `NSC` values as horizontal bars
   ```python
   df = pd.read_csv("out_sensitivity/sensitivity_indices.csv")
   df.plot.barh(x='param_name', y='NSC', title='Sensitivity Indices')
   ```

2. **Prediction Interval Bands**: Plot 70% PI over time for each config
   ```python
   df = pd.read_csv("out_sensitivity/config_pi70.csv")
   for cfg in df['config_id'].unique():
       subset = df[df['config_id'] == cfg]
       plt.fill_between(subset['week'], subset['q15'], subset['q85'], alpha=0.3)
       plt.plot(subset['week'], subset['median'], label=cfg)
   ```

3. **RMSE Comparison**: Bar chart of mean RMSE per configuration

4. **Parameter Sweep**: If extending to multiple delta values, plot NSC vs. delta

## Technical Notes

### Common Random Numbers (CRN)

CRN is enforced by:
1. **Same base_seed across all configs**
2. **Same replicate_id values (0, 1, 2, ...) for corresponding replicates**
3. **Deterministic RNG initialization**: `seed = base_seed + replicate_id`

This ensures that random event streams (e.g., human mobility, mosquito emergence) are synchronized, so differences in outputs are due to parameter changes, not random noise.

### Memory & Performance

- Each run generates ~1 MB of output (weekly_cases + results.csv)
- 15 runs (5 configs × 3 reps) ≈ 15 MB
- 100 runs (5 configs × 20 reps) ≈ 100 MB
- No parallelization yet; runs are sequential

To parallelize:
- Use GNU Parallel or a job scheduler (SLURM, PBS)
- Modify `run_sensitivity.py` to submit jobs instead of running directly

### Extending the Parameter Set

To add a new parameter (e.g., `gamma` - recovery rate):

1. Add to `Params` struct in `MyModel.h`:
   ```cpp
   double gamma{0.14};  // 1/7 day recovery
   ```

2. Use in model code (e.g., `MyHuman.cpp`):
   ```cpp
   if (timeSinceInfection >= (1.0 / params_->gamma)) { ... }
   ```

3. Update `Main.cpp` to read from props:
   ```cpp
   params.gamma = (props.getProperty("gamma") != "") ? 
       repast::strToDouble(props.getProperty("gamma")) : 0.14;
   ```

4. Add to `model.props`:
   ```
   gamma=0.14
   ```

5. Add to `PERTURB_PARAMS` in `run_sensitivity.py`:
   ```python
   PERTURB_PARAMS = ["sigma_M", "sigma_H", "z", "beta_mh", "gamma"]
   ```

## Troubleshooting

### Error: "Cannot open observed data file"

**Cause**: `obs_csv` path is incorrect or file doesn't exist.

**Fix**: Either:
- Provide correct path in `model.props`: `obs_csv=view/incidence_rates_santa_marta.csv`
- Or leave empty to disable RMSE: `obs_csv=`

### Error: "weekly_neighborhood_cases.csv not found"

**Cause**: Model didn't run long enough to generate weekly data.

**Fix**: Check `stop.at` in `model.props`. For 52 weeks, set `stop.at=365`.

### Warning: "Model exited with code 1"

**Cause**: Simulation crashed (e.g., segfault, assertion failure).

**Fix**: Check `simulation.log` in the run directory. Common issues:
- Memory errors (reduce `count.of.humans`)
- Invalid parameter values (negative rates, etc.)

### Sensitivity indices are all near zero

**Possible causes**:
1. **Perturbation too small**: Increase `DELTA` to 0.20 or 0.30
2. **Stochastic noise dominates**: Increase `--replicates` to 10 or 20
3. **Parameter not influential**: Expected for some parameters

## References

- **Saltelli et al. (2008)**: *Global Sensitivity Analysis: The Primer*
- **Common Random Numbers**: Law & Kelton, *Simulation Modeling and Analysis*
- **One-At-a-Time vs. Global SA**: Morris (1991), Sobol (1993)

## Support

For questions or issues:
1. Check this README
2. Review `calibration_I0.py` for similar patterns
3. Inspect `out_sensitivity/meta.json` for run configuration
4. Check model logs in each replicate directory

---

**Last Updated**: 2025-10-29  
**Version**: 1.0  
**Author**: Generated via AI-assisted refactoring

