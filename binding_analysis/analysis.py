import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from statsmodels.api import OLS, add_constant
from statsmodels.stats.diagnostic import acorr_ljungbox, acorr_breusch_godfrey, het_white
from models import model_definitions
from utils import save_combined_csv
from utils import autocorrelation_tests
from utils import validate_data
import traceback
import logging

def process_csv_files_in_folder(input_folder, output_folder):
    for filename in os.listdir(input_folder):
        if not filename.endswith(".csv"):
            continue

        full_path = os.path.join(input_folder, filename)
        logging.info(f"Processing: {full_path}")

        # Use context manager + pandas
        try:
            df = pd.read_csv(full_path)
        except Exception as e:
            logging.error(f"Failed to read {filename}: {e}")
            continue
    
        # Validate required columns and check for nulls
        if not validate_data(df):
            logging.error(f"Skipping {filename} due to validation error.")
            continue

        H0, G0, d_delta_exp = np.genfromtxt(full_path, delimiter=',', names=True, unpack=True)
        H0 = np.array(H0, dtype=float)
        G0 = np.array(G0, dtype=float)
        delta = np.array(d_delta_exp, dtype=float)
        d_delta_exp = np.abs(delta - delta[0])
        d_delta_exp[0] = 0

        models = model_definitions(H0, G0, d_delta_exp)
        results = {}
        output_rows = []

        for model_name, model in models.items():
            try:
                guess = model['initial_guess']
                bounds = model['bounds']
                func = model['lambda']

                params, cov = curve_fit(func, H0, d_delta_exp, p0=guess, bounds=bounds, maxfev=100000)
                fit_vals = func(H0, *params)
                residuals = fit_vals - d_delta_exp

                std_err = np.sqrt(np.diag(cov))
                rss = np.sum(residuals**2)
                r2 = 1 - rss / np.sum((d_delta_exp - np.mean(d_delta_exp))**2)
                n = len(d_delta_exp)
                p = len(params)
                ll = -n / 2 * np.log(rss / n)
                aic = 2 * p - 2 * ll
                bic = p * np.log(n) - 2 * ll
                rmse = np.sqrt(rss / n)

                results[model_name] = (fit_vals, residuals)

                autocorr = autocorrelation_tests(H0, residuals, model_name)

                output_rows.append({
                    "file": filename,
                    "model": model_name,
                    "parameters": params.tolist(),
                    "standard_errors": std_err.tolist(),
                    "r_squared": r2,
                    "AIC": aic,
                    "BIC": bic,
                    "RMSE": rmse,
                    "confidence_intervals": [(v - 1.96*s, v + 1.96*s) for v, s in zip(params, std_err)],
                    "fitted_values": fit_vals.tolist(),
                    "residuals": residuals.tolist(),
                    "H_over_G": (H0 / G0).tolist(),
                    "ljung_stat": autocorr["ljung_stat"],
                    "ljung_p": autocorr["ljung_p"],
                    "ljung_failed": autocorr["ljung_failed"],
                    "bg_test": autocorr["bg_name"],
                    "bg_stat": autocorr["bg_stat"],
                    "bg_p": autocorr["bg_p"],
                    "bg_failed": autocorr["bg_failed"]
                })

                logging.info(f"{model_name} → R²={r2:.3f}, RMSE={rmse:.3f}")

            except Exception as e:
                logging.error(f"Exception in model {model_name}: {e}")
                logging.error(traceback.format_exc())  # <- full traceback
                continue

        output_file = os.path.join(output_folder, filename.replace(".csv", "_results.csv"))
        save_combined_csv(output_rows, output_file)

        plot_path = os.path.join(output_folder, filename.replace(".csv", "_plot.png"))
        plot_results(H0, G0, d_delta_exp, results, plot_path)

def plot_results(H0, G0, d_delta_exp, model_results, filename):
    fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(14, 12))
    model_names = list(model_results.keys())

    for i, model in enumerate(model_names):
        fitted, residuals = model_results[model]
        x = H0 / G0

        axes[i, 0].scatter(x, d_delta_exp, label='Experimental')
        axes[i, 0].plot(x, fitted, color='red', label=f'{model} Fit')
        axes[i, 0].set_title(f'{model} Fit')
        axes[i, 0].legend()

        axes[i, 1].scatter(x, residuals, color='red')
        axes[i, 1].axhline(0, linestyle='--')
        axes[i, 1].set_title(f'{model} Residuals')

    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
