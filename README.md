# BEAST - Binding Evaluation and Analysis Software Tool

This Python-based analysis tool is designed for evaluating titration data from NMR experiments. It fits experimental chemical shift data to multiple thermodynamic binding models and estimates association constants. The tool supports batch analysis, model fitting, statistical diagnostics, and reproducible outputs.

![BEAST](example_data/BEAST_NEW.jpg)

---

## Features

- **Model Support**:
  - Includes `1:1 (HG)`, `1:2 (HG + HG₂)`, `2:1 (HG + H₂G)`, `Dimer (HG + H₂)`, and `Multi-equilibrium (H₂G + HG + H₂)` binding models
  - Model logic is modularized – easily extendable via the `models.py` module

- **Flexible Fitting & Configuration**:
  - Curve fitting using `scipy.optimize.curve_fit` with enhanced initial parameter estimation
  - Smart initial guess refinement strategy to avoid poor local minima
  - All key settings (initial guesses, bounds, residual test thresholds) configurable via `config.yaml`

- **Detailed Fit Diagnostics** (grouped and structured):
  - **Criteria**:
    - Coefficient of determination (R²)
    - Root Mean Square Error (RMSE)
    - Weighted Root Mean Square Error (wRMSE)
    - Akaike Information Criterion (AIC)
    - Bayesian Information Criterion (BIC)
  - **Core Tests**:
    - Shapiro-Wilk (normality test)
    - Ljung-Box test (autocorrelation)
    - Ramsey RESET test (model misspecification)
    - Lag-1 Pearson & Spearman correlation of residuals (PC & SC)
    - Rolling window R² analysis (local structure detection)
  - **Optional Tests**:
    - Run ratio (sign clustering via runs test)
    - Zero-crossings count (ZC)
    - Cook's Distance (outlier/influence detection)
    - Breusch-Godfrey test (higher-order autocorrelation)

- **Visualization**:
  - Fitted curve vs experimental Δδ
  - Raw and normalized residuals
  - Multi-model comparison plots with color coding

- **Structured Output**:
  - Results saved in CSV format (per model, per point, and combined)
  - Diagnostic plots saved as PNG
  - Human-readable diagnostic tables in console and log

- **User Experience**:
  - Logging to terminal and `results/log.txt`
  - Fully modular code structure for reuse and testing
  - CLI interface with customizable options for quick runs
  - Google Colab compatible for web/cloud use

---

## How to Use

### Option 1: Google Colab

Launch in [Google Colab](https://colab.research.google.com/github/Deam0on/B_E_A_S_T/blob/main/example_data/colab_template.ipynb)

A set of experimental data can be found under `example_data/`

### Option 2: Run Locally

1. Clone the repository:

```bash
git clone https://github.com/Deam0on/B_E_A_S_T.git
cd B_E_A_S_T
```

2. Install dependencies:

```bash
pip install -r binding_analysis/requirements.txt
```

3. Place `.csv` files into the `data_input/` directory (or specify a custom directory).

4. Run the tool:

```bash
python binding_analysis/binding_analysis_tool.py [OPTIONS]
```

### Available Options:

| Argument            | Description |
|---------------------|-------------|
| `--config PATH`     | Path to custom `config.yaml` file (default: `binding_analysis/config.yaml`) |
| `--input_dir DIR`   | Override input folder path (default: `data_input`) |
| `--output_dir DIR`  | Override output folder path (default: `results`) |
| `--skip_tests`      | Disable advanced residual diagnostics |
| `--no_normalized`   | Suppress normalized residual plots from output |

### Example:

```bash
python binding_analysis/binding_analysis_tool.py \
  --input_dir example_data \
  --output_dir results \
  --config binding_analysis/config.yaml
```

5. Review the output in the `results/` folder:
   - `*_results.csv`: fit parameters and stats per model
   - `*_combined.csv`: aggregated results across all models
   - `*_plot.png`: fit curves and residuals
   - `log.txt`: execution log

---

## Input File Format

Example data files are provided in the `example_data/` folder.

For detailed input format requirements, see [Input Format and File Structure](https://github.com/Deam0on/B_E_A_S_T/wiki/Input_and_File_Structure) wiki page.

---

## Binding Models

This tool supports multiple binding models for analyzing host-guest interactions:

- **1:1 (HG)** - Simple host-guest binding
- **1:2 (HG + HG₂)** - Sequential binding of two guests
- **2:1 (HG + H₂G)** - Sequential binding of two hosts
- **Dimer (HG + H₂)** - Host-guest binding with host dimerization
- **Multi-equilibrium (H₂G + HG + H₂)** - Complex system with multiple equilibria

For detailed mathematical derivations and model descriptions, see:
[Binding Models and Theory](https://github.com/Deam0on/B_E_A_S_T/wiki/Binding_Models_and_Theory) wiki page

---

## Fit Diagnostics

A comprehensive suite of statistical diagnostics is applied to evaluate model quality:

- **Criteria** - Standard fit quality metrics (R², RMSE, AIC, BIC)
- **Core Tests** - Essential statistical validation tests
- **Optional Tests** - Advanced pattern detection and outlier analysis

See the [Fit Diagnostics](https://github.com/Deam0on/B_E_A_S_T/wiki/Fit_Diagnostics) wiki page for test details, acceptable values, and interpretation.

---

## Output Summary

Each run creates output in the specified results directory:

- **CSV files**: Per-sample results with fit parameters, residuals, metrics, and diagnostic flags
- **Combined CSV**: Aggregated results across all models for easy comparison
- **PNG plots**: Data visualization with fitted curves and residual analysis
- **log.txt**: Comprehensive execution log with model evaluation details

---

## Configuration

The tool's behavior can be customized via `config.yaml`. See the [Configuration Guide](https://github.com/Deam0on/B_E_A_S_T/wiki/Configuration_Guide) wiki page for details on:

- Model-specific parameters and bounds
- Diagnostic test thresholds
- Output preferences
- Advanced fitting options

---

## Additional Resources

- [Project Wiki](https://github.com/Deam0on/B_E_A_S_T/wiki) - Complete documentation
- [CLI Usage and Flags](https://github.com/Deam0on/B_E_A_S_T/wiki/CLI_Usage_and_Flags) - Detailed command-line reference
- [Google Colab Integration](https://github.com/Deam0on/B_E_A_S_T/wiki/Google_Colab_Integration) - Cloud-based usage guide

For questions or contributions, open an issue or pull request on the [GitHub repo](https://github.com/Deam0on/B_E_A_S_T).

---

## Citation

If you use BEAST in your research, please see [CITATION.cff](CITATION.cff) for more details.

---

## License

This project is released under the MIT License. See [LICENSE](LICENSE) for details.

---

## Support this project

You can support development via [Liberapay](https://liberapay.com/Deamoon)
