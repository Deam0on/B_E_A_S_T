
# 🧪 Reproducibility Report

This document provides the necessary information to ensure that the results produced using this repository can be reproduced reliably across different environments.

---

## 📦 Python Environment

This project was developed and tested with:

- Python version: 3.11
- Operating System: Ubuntu 22.04 / Google Colab (Debian-based)

---

## 📚 Dependencies

All required packages are listed in [`requirements.txt`](requirements.txt).  
You can install them using:

```bash
pip install -r requirements.txt
```

To replicate the exact environment used in the original analysis, use the frozen lock file (if present):

```bash
pip install -r requirements_lock.txt
```

To generate a lock file from your current environment:

```bash
pip freeze > requirements_lock.txt
```

---

## 🔁 How to Reproduce the Results

1. Clone the repository:

```bash
git clone https://github.com/Deam0on/mysak_delta_iso.git
cd mysak_delta_iso
```

2. Create a virtual environment (recommended):

```bash
python -m venv env
source env/bin/activate  # On Windows: env\Scripts\activate
```

3. Install dependencies:

```bash
pip install -r binding_analysis/requirements.txt
```

4. Run the script:

```bash
cd binding_analysis
python binding_analysis_tool.py
```

5. Upload your CSV file when prompted or use `example_data/test.csv` to test functionality.

---

## 🧪 Sample Dataset

An example input file is included in:

```bash
example_data/test.csv
```

You can copy or move it into the `data_input/` directory before running the tool.

---

## 📝 Logging and Traceability

- All logs and outputs will be saved in the `results/` folder.
- Logs include model parameters, metrics, autocorrelation diagnostics, and error messages if any.

---

## 🧠 Notes

- The models used involve numerical optimization and may produce slightly different fits depending on environment (especially if run on CPU vs GPU or different BLAS libraries).
- Ensure reproducibility by fixing the initial guesses or setting seeds where randomness is introduced (e.g., in optimization).

