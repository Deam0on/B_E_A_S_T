
# Binding Isotherm Analysis Tool

This Python-based analysis tool is designed for evaluating titration data from NMR experiments. It fits experimental chemical shift data to multiple thermodynamic binding models and estimates association constants. The tool supports batch analysis, model fitting, statistical diagnostics, and reproducible outputs.

---

## Features

- Supports multiple host-guest binding models:
  - 1:1, 1:2, 2:1, dimerization, and multi-equilibrium
- Configurable fitting parameters via `config.yaml`
- Curve fitting using `scipy.optimize.curve_fit`
- Initial guess optimization and bound control
- Fit diagnostics including:
  - R², AIC, BIC, RMSE
  - Ljung-Box and Breusch-Godfrey or White's test for residual autocorrelation
- Automatic plot generation (fit and residuals)
- Structured output in CSV and PNG format
- Complete logging to terminal and `results/log.txt`
- Example dataset and reproducibility instructions included
- Google Colab support for cloud-based execution

---

## Input File Format

Input `.csv` files must be placed in the `data_input/` folder. Each file must contain the following columns:

| Column  | Description                        |
|---------|------------------------------------|
| `H`     | Host concentration (mol/L)         |
| `G`     | Guest concentration (mol/L)        |
| `delta` | Observed chemical shift (Hz)       |

Do not rename these column headers.

An example input file is provided in `example_data/test.csv`.

---

## How to Use

### Option 1: Run in Google Colab

Launch an interactive session:

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Deam0on/mysak_delta_iso/blob/main/example_data/colab_template.ipynb)

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

3. Place your `.csv` files into the `data_input/` folder.

4. Run the tool:

```bash
python binding_analysis_tool.py
```

5. Review the output in the `results/` folder:
   - Fit statistics and parameters in `*_results.csv`
   - Fit and residual plots in `*_plot.png`
   - Full execution log in `log.txt`

---

## Configuration

The `config.yaml` file allows full control of:

- Initial parameter guesses
- Parameter bounds
- Residual test lags
- Max iterations for curve fitting
- Input/output directories

Modify `binding_analysis/config.yaml` to suit your experiment.

---

## Output Summary

Each model, for each dataset, produces the following metrics:

- Optimized parameters with standard errors
- 95% confidence intervals
- R² (coefficient of determination)
- AIC and BIC model selection criteria
- RMSE (root mean squared error)
- Residual diagnostics:
  - Ljung-Box test (autocorrelation)
  - Breusch-Godfrey or White's test (fallback)

All metrics are written into a structured CSV file.

---

## Folder Structure

```
mysak_delta_iso/
├── binding_analysis/
│   ├── binding_analysis_tool.py
│   ├── analysis.py
│   ├── models.py
│   ├── utils.py
│   ├── config.yaml
│   └── requirements.txt
├── example_data/
│   ├── colab_template.ipynb
│   └── test.csv
├── reproducibility.md
└── README.md
```

---

## Supported Models

| Model  | Description                       |
|--------|-----------------------------------|
| 1:1    | Single host-guest association     |
| 1:2    | Host binds two guest molecules    |
| 2:1    | Two hosts bind one guest          |
| Dimer  | Dimerization process              |
| Multi  | Multi-equilibrium with H, G, HG, H₂G species |

---

## Reproducibility

Refer to `reproducibility.md` for full environment information and step-by-step instructions to ensure consistent results across systems.

---

## Dependencies

All dependencies are listed in `requirements.txt`:

- numpy
- pandas
- matplotlib
- scipy
- statsmodels
- sympy
- pyyaml

Install using:

```bash
pip install -r requirements.txt
```

---

## License

This project is released under the MIT License. See `LICENSE` for details.
