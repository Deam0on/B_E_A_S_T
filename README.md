# Binding Isotherm Analysis Tool

This Python-based analysis tool is designed for evaluating titration data from NMR experiments. It fits experimental chemical shift data to multiple thermodynamic binding models and estimates association constants. The tool supports batch analysis, model fitting, statistical diagnostics, and reproducible outputs.

---

## Features

- Supports multiple host-guest binding models:
  - 1:1, 1:2, 2:1, dimerization, and multi-equilibrium models
- Full decoupling of model logic from analysis code
- Curve fitting using `scipy.optimize.curve_fit`
- Smart initial guess optimizer (gradient-based local search)
- Automatic model selection based on AIC, BIC, and RMSE
- Statistical diagnostics:
  - R², AIC, BIC, RMSE
  - Ljung-Box and Breusch-Godfrey/White's test for residual autocorrelation
- Plot generation (fit + residuals)
- Structured outputs in CSV and PNG format
- CLI mode with automatic model comparison
- Unit test suite for all models using `pytest`
- Google Colab support for cloud-based execution

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
## Initial Guess Handling

Each model in the tool defines its own **initial guess** and **parameter bounds** directly within the `model_definitions(...)` structure.

### Per-Model Configuration Example

```python
"1:2": {
    "function": binding_isotherm_1_2,
    "initial_guess": [100, 100, 100, 100],
    "bounds": ([0, 0, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf]),
    ...
}
```

These defaults ensure that the fitting begins from a reasonable parameter space.

### Smart Optimization Logic

Before fitting, the code automatically performs **initial guess refinement** using an iterative optimizer:

```python
guess = optimize_initial_guess(
    model["lambda"],
    model["initial_guess"],
    model["bounds"],
    H0,
    G0,
    d_delta_exp,
)
```

This function:

- Starts from the provided guess
- Tweaks each parameter up/down (step = 10 by default)
- Accepts the change only if it improves the RMSE (Root Mean Square Error)
- Stops early if no further improvement is found

This improves convergence while still respecting model-specific parameter structures.

### Benefits

- Robust against poor initial guesses
- Reduces manual tuning of `initial_guess`
- Keeps fitting logic modular and reusable

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
