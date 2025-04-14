# Binding Isotherm Analysis Tool

This Python-based analysis tool is designed for evaluating titration data from NMR experiments. It fits experimental chemical shift data to multiple thermodynamic binding models and estimates association constants. The tool supports batch analysis, model fitting, statistical diagnostics, and reproducible outputs.

---

## Features (Section under construction)

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
[Open in Colab](https://colab.research.google.com/github/Deam0on/mysak_delta_iso/blob/custom_beta/example_data/colab_template.ipynb)

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

---

## Input File Format

See [Input Format and File Structure](Input_and_File_Structure)

---

## Binding Models — Equations and Definitions

This repository includes multiple binding models used to fit NMR titration data. Below are the corresponding equations and definitions of parameters.

For detailed model description, see [Binding Models and Theory](Binding_Models_and_Theory)

---

## Residual Diagnostics

A rich suite of statistical diagnostics is applied to each model's residuals, grouped into:

- **Criteria metrics** (e.g., R², AIC)
- **Main diagnostics** (e.g., Ljung-Box, RESET)
- **Optional residual pattern checks** (e.g., Pearson/Spearman correlation, spectral energy ratio)

These are rendered into a per-model table for easy comparison.

See the [Residual Diagnostics wiki page](https://github.com/Deam0on/mysak_delta_iso/wiki/Residual-Diagnostics) for test details, acceptable values, and interpretation.

---

## Output Summary

Each run creates a `results/` folder with:

- Per-sample **CSV summaries** (with fit, residuals, metrics, flags)
- Per-sample **plots**: data, residuals, and optionally normalized residuals
- A comprehensive **log.txt** with step-by-step model evaluation

The [Output Format wiki page](https://github.com/Deam0on/mysak_delta_iso/wiki/Output-Files) describes all output files and CSV structure.

---
## Additional Resources

- 📚 [Project Wiki](https://github.com/Deam0on/mysak_delta_iso/wiki)
- 📂 [Example Dataset](https://github.com/Deam0on/mysak_delta_iso/tree/main/data_input)
- 📈 [Adding New Models](https://github.com/Deam0on/mysak_delta_iso/wiki/Adding-New-Models)

For questions or contributions, open an issue or pull request on the [GitHub repo](https://github.com/Deam0on/mysak_delta_iso).

---

## License

This project is released under the MIT License. See `LICENSE` for details.

## Support this project

You can support me via [liberpay](liberapay.com/Deamoon)
