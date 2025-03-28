
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


## Input File Format

Input `.csv` files must be placed in the `data_input/` folder and contain the following columns:

| Column | Description                      |
|--------|----------------------------------|
| `H`    | Host concentration               |
| `G`    | Guest concentration              |
| `delta`| Observed chemical shift (Hz)     |

> Do not rename these columns – they are required for parsing.

Example dataset is included in `example_dataset/test.csv`

---

## How to Run

You can run a template:
> [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Deam0on/mysak_delta_iso/blob/main/example_data/colab_template.ipynb)



Or manually:
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
## Folder Structure

```
mysak_delta_iso/
├── binding_analysis/
│   ├── binding_analysis_tool.py
│   ├── analysis.py
│   ├── utils.py
│   ├── models.py
│   └── requirements.txt
├── colab_template.ipynb
└── README.md
```

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
