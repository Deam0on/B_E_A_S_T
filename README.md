
# Binding Isotherm Analysis Tool

This Python tool analyzes titration experiment data to determine association constants using multiple isotherm models. It fits experimental chemical shift data to binding models such as 1:1, 1:2, 2:1, dimer, and multi-binding equilibria.

---

## Features

- Multiple binding models
- Automatic parameter optimization
- Curve fitting using `scipy.optimize.curve_fit`
- Fit statistics: R², AIC, BIC, RMSE
- Residual analysis with Ljung-Box and Breusch-Godfrey/White’s test
- Fit and residual plots saved per dataset
- Batch processing of CSV files
- Logging of progress and results

---

## Folder Structure

```
binding_analysis/
├── binding_analysis_tool.py       # Main script (entry point)
├── models.py                      # Model definitions and initial guess logic
├── analysis.py                    # Curve fitting, evaluation, and plotting
├── utils.py                       # CSV handling, cleanup, logging setup
├── requirements.txt               # All required dependencies
├── README.md                      # This file
├── data_input/                    # Folder for input CSV files
└── results/                       # Output plots and results
```

---

## Input File Format

Input `.csv` files must be placed in the `data_input/` folder and contain the following columns:

| Column | Description                      |
|--------|----------------------------------|
| `H`    | Host concentration               |
| `G`    | Guest concentration              |
| `delta`| Observed chemical shift (Hz)     |

> Do not rename these columns – they are required for parsing.

---

## How to Run

1. Install the dependencies:

```bash
pip install -r requirements.txt
```

2. Add your `.csv` data files into the `data_input/` directory.

3. Run the tool:

```bash
python binding_analysis_tool.py
```

4. Output will be saved to the `results/` folder:
   - Fitted values and residuals (`*_results.csv`)
   - Fit and residual plots (`*_plot.png`)
   - Logging output printed to console

---

## Supported Binding Models

| Model  | Description                   |
|--------|-------------------------------|
| 1:1    | Simple 1:1 host-guest         |
| 1:2    | Host to two guest molecules   |
| 2:1    | Two hosts to one guest        |
| Dimer  | Dimerizing system             |
| Multi  | Multi-equilibrium (H, G, HG, H₂G) |

---

## Output Summary (per model per file)

- Optimized parameters
- Standard errors
- 95% confidence intervals
- R², AIC, BIC, RMSE
- Residual diagnostics (autocorrelation tests)

---

## Dependencies

- numpy
- pandas
- matplotlib
- scipy
- statsmodels
- sympy

Install them using:

```bash
pip install -r requirements.txt
```
