# Binding Isotherm Analysis Tool

This Python-based analysis tool is designed for evaluating titration data from NMR experiments. It fits experimental chemical shift data to multiple thermodynamic binding models and estimates association constants. The tool supports batch analysis, model fitting, statistical diagnostics, and reproducible outputs.

---

## Features

- Supports multiple host-guest binding models:
  - 1:1, 1:2, 2:1, dimerization, and multi-equilibrium
- Curve fitting using `scipy.optimize.curve_fit`
- Configurable fitting parameters via `config.yaml`
- Fit diagnostics including:
  - R², AIC, BIC, RMSE
  - Ljung-Box and Breusch-Godfrey or White's test for residual autocorrelation
- Plot generation (fit and residuals)
- Structured output in CSV and PNG format
- Logging to terminal and `results/log.txt`
- Reproducibility guidance in `REPRODUCIBILITY.md`
- Unit test suite for all models using pytest
- Google Colab support for cloud-based execution

---

## Input File Format

Input `.csv` files must be placed in the `data_input/` folder. Each file must contain the following columns:

| Column  | Description                        |
|---------|------------------------------------|
| `H`     | Host concentration (mol/L)         |
| `G`     | Guest concentration (mol/L)        |
| `delta` | Observed chemical shift (Hz)       |

An example input file is provided in `example_data/test.csv`.

---

## Supported Binding Models

### 1:1 Binding

**Equilibrium:**

H + G ⇌ HG  
Ka = [HG] / ([H][G])

**Observed shift:**

δ_obs = δ_G + (δ_HG - δ_G) × [HG] / [G]_0

---

### 1:2 Binding

**Equilibria:**

H + G ⇌ HG (K₁)  
HG + G ⇌ HG₂ (K₂)

**Observed shift:**

δ_obs = (δ_G × [G]_0 + δ_HG × [HG] + δ_HG₂ × [HG₂]) / [G]_0

---

### 2:1 Binding

**Equilibria:**

H + G ⇌ HG (K₁)  
H + HG ⇌ H₂G (K₂)

**Observed shift:**

δ_obs = (δ_G × [G]_0 + δ_HG × [HG] + δ_H₂G × [H₂G]) / [G]_0

---

### Dimerization

**Equilibria:**

H + G ⇌ HG  
2H ⇌ H₂ (K_d)

**Observed shift:**

δ_obs = (δ_G / (1 + K_a[H]) + δ_HG × K_a[H]) / (1 / (1 + K_a[H]) + K_a[H])

---

### Multi-binding (HG and H₂G)

**Equilibria:**

H + G ⇌ HG (K_HG)  
2H ⇌ H₂ (K_d)  
H₂ + G ⇌ H₂G (K_H₂G)

**Observed shift:**

δ_obs = (δ_G × [G] + δ_HG × [HG] + δ_H₂G × [H₂G]) / [G]_0

---

## How to Use

### Option 1: Google Colab

Launch in Colab:  
[Open in Colab](https://colab.research.google.com/github/Deam0on/mysak_delta_iso/blob/main/example_data/colab_template.ipynb)

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
python binding_analysis_tool.py
```

5. Review the output in the `results/` folder:
   - `*_results.csv`: fit parameters and stats
   - `*_plot.png`: fit and residuals
   - `log.txt`: execution log

---

## Configuration

Edit `binding_analysis/config.yaml` to adjust:

- Initial parameter guesses
- Parameter bounds
- Input/output folder paths
- Residual test lag count
- Fitting settings (e.g., max iterations)

---

## Testing

A full test suite is included in `tests/test_models.py`.

To run locally:

```bash
pytest tests/
```

Tests ensure:

- Model functions return finite values
- curve_fit optimization converges
- Estimated parameters are within tolerance of known values

Tests are also executed automatically via GitHub Actions.

To check test coverage:

```bash
pytest --cov=binding_analysis --cov-report=term
```

---

## Output Summary

Each model-dataset pair produces:

- Optimized parameters with standard errors
- 95% confidence intervals
- R², AIC, BIC, RMSE
- Residual diagnostics
- Output CSV and PNG plots
- Combined execution log

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
├── results/
├── tests/
│   └── test_models.py
├── REPRODUCIBILITY.md
└── README.md
```

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

## Reproducibility

All runtime and version info for reproducibility is in `REPRODUCIBILITY.md`.

---

## License

This project is released under the MIT License. See `LICENSE` for details.