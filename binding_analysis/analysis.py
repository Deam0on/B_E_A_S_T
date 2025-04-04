import os
import numpy as np
import pandas as pd
import traceback
import logging
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from statsmodels.api import OLS, add_constant
from statsmodels.stats.diagnostic import acorr_ljungbox, acorr_breusch_godfrey, het_white

from models import model_definitions
from utils import save_combined_csv, autocorrelation_tests, validate_data

def smart_initial_guess(model_func, guess, bounds, H0, d_delta_exp, step=10, max_iter=10):
    best_guess = np.array(guess)
    best_rmse = np.inf

    for _ in range(max_iter):
        improved = False
        for i in range(len(best_guess)):
            for delta in [-step, step]:
                candidate = best_guess.copy()
                candidate[i] += delta
                candidate[i] = np.clip(candidate[i], bounds[0][i], bounds[1][i])
                try:
                    fit_vals = model_func(H0, *candidate)
                    rmse = np.sqrt(np.mean((fit_vals - d_delta_exp) ** 2))
                    if rmse < best_rmse:
                        best_guess = candidate
                        best_rmse = rmse
                        improved = True
                except:
                    continue
        if not improved:
            break

    return best_guess

def compare_models_by_metric(output_rows, metric="AIC"):
    sorted_models = sorted(output_rows, key=lambda r: r[metric])
    logging.info("\nModel ranking by %s:" % metric)
    for rank, row in enumerate(sorted_models, 1):
        logging.info(f"{rank}. {row['model']} | {metric} = {row[metric]:.2f} | RMSE = {row['RMSE']:.4f}")
    return sorted_models

def advanced_residual_diagnostics(H0, residuals, model_name):
    autocorr = autocorrelation_tests(H0, residuals, model_name)

    residuals_std = (residuals - np.mean(residuals)) / np.std(residuals)
    normality_pass = np.all(np.abs(residuals_std) < 3)

    skewness = pd.Series(residuals).skew()
    kurtosis = pd.Series(residuals).kurtosis()

    logging.info(f"Additional diagnostics for {model_name}:")
    logging.info(f"Skewness: {skewness:.3f}, Kurtosis: {kurtosis:.3f}")
    if not normality_pass:
        logging.warning("Residuals may not be normally distributed (outliers or heavy tails).")

    return {
        **autocorr,
        "skewness": skewness,
        "kurtosis": kurtosis,
        "normality_pass": normality_pass
    }

def process_csv_files_in_folder(config):
    input_folder = config["general"]["input_dir"]
    output_folder = config["general"]["results_dir"]
    maxfev = config["general"]["maxfev"]
    lags = config["general"]["lags"]

    for filename in os.listdir(input_folder):
        if not filename.endswith(".csv"):
            continue

        full_path = os.path.join(input_folder, filename)
        logging.info(f"Processing: {full_path}")

        try:
            df = pd.read_csv(full_path)
        except Exception as e:
            logging.error(f"Failed to read {filename}: {e}")
            continue

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
                logging.info(f"Evaluating model: {model_name}")
                guess = smart_initial_guess(model['lambda'], model['initial_guess'], model['bounds'], H0, d_delta_exp)
                bounds = model['bounds']
                func = model['lambda']

                params, cov = curve_fit(func, H0, d_delta_exp, p0=guess, bounds=bounds, maxfev=maxfev)
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

                diagnostics = advanced_residual_diagnostics(H0, residuals, model_name)

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
                    **diagnostics
                })

                logging.info(f"Model {model_name} fit completed")
                logging.info(f"Parameters: {params}")
                logging.info(f"R²: {r2:.4f}, AIC: {aic:.2f}, BIC: {bic:.2f}, RMSE: {rmse:.4f}")

            except Exception as e:
                logging.error(f"Exception in model {model_name}: {e}")
                logging.error(traceback.format_exc())
                continue

        sorted_models = compare_models_by_metric(output_rows, metric="AIC")

        output_file = os.path.join(output_folder, filename.replace(".csv", "_results.csv"))
        save_combined_csv(sorted_models, output_file)

        plot_path = os.path.join(output_folder, filename.replace(".csv", "_plot.png"))
        plot_results(H0, G0, d_delta_exp, results, plot_path)

def plot_results(H0, G0, d_delta_exp, model_results, filename):
    fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(14, 12))
    model_names = list(model_results.keys())

    for i, model in enumerate(model_names):
        fitted, residuals = model_results[model]
        x = H0 / G0

        # Left plot: Fit
        axes[i, 0].scatter(x, d_delta_exp, label='Experimental')
        axes[i, 0].plot(x, fitted, color='red', label=f'{model} Fit')
        axes[i, 0].set_title(f'{model} Fit')
        axes[i, 0].set_xlabel('[H]/[G]')
        axes[i, 0].set_ylabel('Δδ [Hz]')
        axes[i, 0].legend()

        # Right plot: Residuals
        axes[i, 1].scatter(x, residuals, color='red')
        axes[i, 1].axhline(0, linestyle='--')
        axes[i, 1].set_title(f'{model} Residuals')
        axes[i, 1].set_xlabel('[H]/[G]')
        axes[i, 1].set_ylabel('Residual [Hz]')

    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
