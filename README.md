# Binding Isotherm Analysis Tool

This Python-based analysis tool is designed for evaluating titration data from NMR experiments. It fits experimental chemical shift data to multiple thermodynamic binding models and estimates association constants. The tool supports batch analysis, model fitting, statistical diagnostics, and reproducible outputs.

---

## Features

- Supports multiple host-guest binding models:
  - 1:1, 1:2, 2:1, dimerization, and multi-equilibrium models
- Full decoupling of model logic from analysis code
- Curve fitting using `scipy.optimize.curve_fit`
- Smart initial guess refinement strategy
- Configurable fitting parameters via `config.yaml`
- Fit diagnostics including:
  - R², AIC, BIC, RMSE, Weighted RMSE
  - Ljung-Box, Breusch-Godfrey or White's test for residual autocorrelation
  - Ramsey RESET Test (for model misspecification)
  - Cook's Distance (influence detection)
  - Skewness, Kurtosis, and D’Agostino’s normality test
- Plot generation:
  - Fitted curves
  - Raw residuals
  - Normalized residuals (Δδ-relative)
- Structured output in CSV and PNG format
- Logging to terminal and `results/log.txt`
- Reproducibility guidance in `REPRODUCIBILITY.md`
- Unit test suite for all models using pytest
- Google Colab support for cloud-based execution

---

## How to Use

### Option 1: Google Colab

Launch in Colab:  
[Open in Colab](https://colab.research.google.com/github/Deam0on/mysak_delta_iso/blob/custom_statistics/example_data/colab_template.ipynb)

### Option 2: Run Locally

1. Clone the repository:

```bash
git clone https://github.com/Deam0on/mysak_delta_iso.git
cd mysak_delta_iso/binding_analysis
```

2. Install dependencies:

```bash
pip install -r requirements.txt
```

3. Place `.csv` files into the `data_input/` directory.

4. Run the tool:

```bash
python binding_analysis_tool.py [OPTIONS]
```

### Available Options:

| Argument            | Description |
|---------------------|-------------|
| `--config PATH`     | Path to custom `config.yaml` file (default: `config.yaml`) |
| `--input_dir DIR`   | Override input folder path (default: `data_input`) |
| `--output_dir DIR`  | Override output folder path (default: `results`) |
| `--skip_tests`      | Disable residual diagnostics like Ljung-Box, BG/White, RESET, Cook’s Distance |
| `--no_normalized`   | Suppress normalized residual plots (Δδ-relative) from output |

### Example:

```bash
python binding_analysis_tool.py \
  --input_dir custom_data \
  --output_dir custom_results \
  --skip_tests \
  --no_normalized
```

5. Review the output in the `results/` folder:
   - `*_results.csv`: fit parameters and stats
   - `*_plot.png`: fit and residuals
   - `log.txt`: execution log

6. Optionally run with automatic model comparison:

```bash
python binding_analysis_tool.py --auto-select
```
---

## Input File Format

See [Input Format and File Structure](Input_and_File_Structure)

---

## Binding Models — Equations and Definitions

This repository includes multiple binding models used to fit NMR titration data. Below are the corresponding equations and definitions of parameters.

For detailed model description, see [Binding Models and Theory](Binding_Models_and_Theory)

---

## Results Interpretation Guide

### Metrics and Diagnostics

| Metric       | Description |
|--------------|-------------|
| **R²**       | Coefficient of determination — closer to 1 is better |
| **AIC/BIC**  | Model selection criteria — lower is better |
| **RMSE**     | Root Mean Squared Error — unnormalized residual size |
| **Weighted RMSE** | RMSE scaled by the Δδ range — good for comparing across datasets |
| **Ljung-Box**| Detects autocorrelation — p < 0.05 is problematic |
| **Breusch-Godfrey / White** | Detects serial correlation / heteroscedasticity |
| **Ramsey RESET** | Checks for non-linearity or omitted variables |
| **Cook’s Distance** | Flags highly influential data points |
| **Skew / Kurtosis / Normality p** | Indicates Gaussian behavior of residuals |

### Cookbook for Model Evaluation

1. **Start with AIC/BIC**: Prefer models with the lowest values.
2. **Verify fit quality using RMSE & Weighted RMSE**: Smaller is better.
3. **Check residual plots**:
   - Look for randomness (no trends)
   - Normalize view via “Normalized Residuals” for cross-experiment comparison
4. **Review Ljung-Box and BG/White Tests**:
   - If both fail → model likely poorly specified
5. **Ramsey RESET test fails?**
   - Try alternate model or include missing terms
6. **Large Cook’s Distance points?**
   - Validate or consider outlier exclusion
7. **Skewness/Kurtosis too high or p < 0.05 in normality test?**
   - Might indicate non-random structure → residuals not white noise

Use these steps as a checklist to pick the best-fitting model per dataset.

---
## Zero-Crossing Similarity

The tool computes the similarity between model residuals and white noise by comparing zero-crossings:

- **Zero-crossings**: Number of sign changes in residuals.
- **Similarity to white noise**: % overlap with simulated normal noise.

This metric is logged and saved in the results file. High similarity (>80%) suggests good model randomness.

---
## Configuration

See [Configuration Guide](Configuration_Guide)

---

## Output Summary

Each model-dataset pair produces:

- Optimized parameters with standard errors
- 95% confidence intervals
- R², AIC, BIC, RMSE
- Residual diagnostics
- Output CSV and PNG plots
- Combined execution log
- Ljung-Box statistics and pass/fail
- BG/White statistics and pass/fail
- Cook’s Distance max and extreme counts
- Ramsey RESET test (if applicable)
- Zero-crossing similarity (%)

---

## Dependencies

Install all using:

```bash
pip install -r requirements.txt
```

- numpy
- pandas
- matplotlib
- scipy
- statsmodels
- sympy
- pyyaml
- pytest
- pytest-cov

---

## License

This project is released under the MIT License. See `LICENSE` for details.

## Support this project

You can support me via [liberpay](liberapay.com/Deamoon)
