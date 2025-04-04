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
[Open in Colab](https://colab.research.google.com/github/Deam0on/mysak_delta_iso/blob/beta/example_data/colab_template.ipynb)

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

Input `.csv` files must be placed in the `data_input/` folder. Each file must contain the following columns:

| Column  | Description                        |
|---------|------------------------------------|
| `H`     | Host concentration (mol/L)         |
| `G`     | Guest concentration (mol/L)        |
| `delta` | Observed chemical shift (Hz)       |

An example input file is provided in `example_data/test.csv`.

---

## Binding Models — Equations and Definitions

This repository includes multiple binding models used to fit NMR titration data. Below are the corresponding equations and definitions of parameters.

### General Notation

- $H_0$: Total concentration of host  
- $G_0$: Total concentration of guest  
- $H_\text{free}, G_\text{free}$: Free (unbound) host and guest concentrations  
- $K_a$: Association constant of first binding event  
- $K_d$: Association constant of second binding event  
- $K_{HG}, K_{H_2G}$: Stepwise association constants  
- $\delta_\text{obs}$: Observed chemical shift  
- $\delta_\text{free}$: Chemical shift of free species  
- $\delta_\infty$: Limiting shift(s) of bound species  
- $\Delta \delta = \delta_\text{obs} - \delta_\text{obs}(H_0 = 0)$: Normalized observed shift  

---

### 1:1 Binding Model

#### Equilibrium

$$
H + G \leftrightarrow HG \quad K_a = \frac{[HG]}{[H][G]}
$$

#### Complex concentration

$$
[HG] = \frac{1}{2} \left( H_0 + G_0 + \frac{1}{K_a} - \sqrt{(H_0 + G_0 + \frac{1}{K_a})^2 - 4 H_0 G_0} \right)
$$

#### Observed shift

$$
\delta_\text{obs} = \delta_\text{free} \cdot \frac{[HG] - G_0}{G_0} + \delta_\infty \cdot \frac{[HG]}{G_0}
$$

---

### 1:2 Binding Model (Host:Guest)

#### Equilibria

$$
\begin{align*}
H + G \leftrightarrow HG \quad K_a = \frac{[HG]}{[H][G]} \\
HG + G \leftrightarrow HG_2 \quad K_d = \frac{[HG_2]}{[HG][G]}
\end{align*}
$$

#### Cubic equation for $G_\text{free}$

$$
aG^3 + bG^2 + cG + d = 0
$$

with:

$$
\begin{align*}
a &= K_a K_d \\
b &= K_a (2K_d H_0 - K_d G_0 + 1) \\
c &= K_a(H_0 - G_0) + 1 \\
d &= -G_0
\end{align*}
$$

#### Observed shift

$$
\delta_\text{obs} = \frac{
\delta_{\infty,1} \cdot G_0 K_a G_\text{free} + 
\delta_{\infty,2} \cdot G_0 K_a K_d G_\text{free}^2
}{
1 + K_a G_\text{free} + K_a K_d G_\text{free}^2
}
$$

---

### 2:1 Binding Model (Host:Guest)

#### Equilibria

$$
\begin{align*}
H + G \leftrightarrow HG \quad K_a = \frac{[HG]}{[H][G]} \\
HG + H \leftrightarrow H_2G \quad K_d = \frac{[H_2G]}{[HG][H]}
\end{align*}
$$

#### Cubic equation for $H_\text{free}$

$$
aH^3 + bH^2 + cH + d = 0
$$

with:

$$
\begin{align*}
a &= K_a K_d \\
b &= K_a (2K_d G_0 - K_d H_0 + 1) \\
c &= K_a(G_0 - H_0) + 1 \\
d &= -H_0
\end{align*}
$$

#### Observed shift

$$
\delta_\text{obs} = \frac{
\delta_{\infty,1} \cdot G_0 K_a H_\text{free} + 
2 \cdot \delta_{\infty,2} \cdot G_0 K_a K_d H_\text{free}^2
}{
G_0 \cdot \left(1 + K_a H_\text{free} + K_a K_d H_\text{free}^2\right)
}
$$

---

### Dimerization Model

#### Equilibria

$$
\begin{align*}
2H \leftrightarrow H_2 \quad K_d = \frac{[H_2]}{[H]^2} \\
H + G \leftrightarrow HG \quad K_a = \frac{[HG]}{[H][G]}
\end{align*}
$$

#### Quadratic equation for $H_\text{free}$

$$
-2K_d H_\text{free}^2 - (K_a G_0 + 1) H_\text{free} + H_0 = 0
$$

#### Observed shift

$$
\delta_\text{obs} = \frac{
\delta_{\infty,1} \cdot \frac{G_0}{1 + K_a H_\text{free}} +
\delta_{\infty,2} \cdot K_a H_\text{free} G_0
}{
\frac{G_0}{1 + K_a H_\text{free}} +
K_a H_\text{free} G_0
}
$$

---

### Multi-Component Model (Dimer, HG, H₂G)

#### Equilibria

$$
\begin{align*}
2H \leftrightarrow H_2 \quad K_d = \frac{[H_2]}{[H]^2} \\
H + G \leftrightarrow HG \quad K_{HG} = \frac{[HG]}{[H][G]} \\
H_2 + G \leftrightarrow H_2G \quad K_{H_2G} = \frac{[H_2G]}{[H_2][G]}
\end{align*}
$$

#### Iterative solution

Solve iteratively for $H_\text{free}$ and $G_\text{free}$ using:

$$
H_0 = H_\text{free} + 2 K_d H_\text{free}^2 + K_{HG} H_\text{free} G_\text{free} + 2 K_{H_2G} K_d H_\text{free}^2 G_\text{free}
$$

$$
G_0 = G_\text{free} + K_{HG} H_\text{free} G_\text{free} + K_{H_2G} K_d H_\text{free}^2 G_\text{free}
$$

#### Observed shift

$$
\delta_\text{obs} = \frac{
\delta_G \cdot K_d H_\text{free}^2 +
\delta_{HG} \cdot K_{HG} H_\text{free} G_\text{free} +
\delta_{H_2G} \cdot K_{H_2G} K_d H_\text{free}^2 G_\text{free}
}{
K_d H_\text{free}^2 +
K_{HG} H_\text{free} G_\text{free} +
K_{H_2G} K_d H_\text{free}^2 G_\text{free}
}
$$

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

Edit `binding_analysis/config.yaml` to adjust:

- Initial parameter guesses
- Parameter bounds
- Input/output folder paths
- Residual test lag count
- Fitting settings (e.g., max iterations)

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
