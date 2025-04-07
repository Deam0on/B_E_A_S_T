import os
import numpy as np
import pandas as pd
import traceback
import logging
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from statsmodels.api import OLS, add_constant
from statsmodels.stats.diagnostic import acorr_ljungbox, acorr_breusch_godfrey, het_white
from utils import collect_global_max_deltadelta
from utils import autocorrelation_tests, validate_data, custom_residual_pattern_test
from scipy.stats import pearsonr, spearmanr
from statsmodels.sandbox.stats.runs import runstest_1samp
from scipy.fft import rfft
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
    from tabulate import tabulate

    sorted_models = sorted(output_rows, key=lambda r: r[metric])
    logging.info("\nModel ranking by %s (lower is better):", metric)

    for rank, row in enumerate(sorted_models, 1):
        model = row["model"]
        aic = row["AIC"]
        bic = row["BIC"]
        r2 = row["r_squared"]
        rmse = row["RMSE"]
        wrmse = row.get("weighted_RMSE", None)
        norm_raw = row.get("normality_pass", True)

        table_data = [
            ["Metric", "Value", "Acceptable Range"],
            ["R²", f"{r2:.4f}", "> 0.98"],
            ["RMSE", f"{rmse:.4f}", "As low as possible"],
            ["Weighted RMSE", f"{wrmse:.4f}" if wrmse is not None else "n/a", "As low as possible"],
            ["AIC", f"{aic:.2f}", "As low as possible"],
            ["BIC", f"{bic:.2f}", "As low as possible"],
            ["Normality test", "✓" if norm_raw else "⚠️", "True / False"],
        ]

        ljung = row.get("ljung_stat")
        if ljung is not None:
            table_data.append(["Ljung-Box stat", f"{ljung:.3f}", "p > 0.05"])

        reset_p = row.get("reset_p")
        if reset_p is not None:
            table_data.append(["RESET p", f"{reset_p:.4f}", "p > 0.05"])

        comp_stats = row.get("composite_stats", {})
        if comp_stats:
            table_data.extend([
                ["Pearson corr", f"{comp_stats['pearson_corr']:.2f}", "|r| < 0.35"],
                ["Spearman corr", f"{comp_stats['spearman_corr']:.2f}", "|r| < 0.35"],
                ["Rolling R²", f"{comp_stats['avg_rolling_r2']:.2f}", "< 0.35"]
                # ["Run ratio", f"{comp_stats['run_ratio']:.2f}", "0.7–1.3"]
            ])

        logging.info(f"\n{rank}. Model: {model}\n" + tabulate(table_data, headers="firstrow", tablefmt="fancy_grid"))

        custom_corr_flagged = row.get("composite_flagged")
        custom_corr_symbol = "✓" if custom_corr_flagged is False else ("⚠️" if custom_corr_flagged is True else "–")

        # if "normality_pass" not in row:
        #     logging.warning(f"normality_pass not found in diagnostics for model {model}")
        # else:
        #     logging.warning(f"Alles OK (normality")
    
        # logging.info(f"{rank}. {model}")
        # logging.info(f"    R² = {r2:.4f} | RMSE = {rmse:.4f}" + (f" | wRMSE = {wrmse:.4f}" if wrmse is not None else ""))
        # logging.info(f"    AIC = {aic:.2f} | BIC = {bic:.2f}")
        # logging.info(f"    Skewness = {row.get('skewness', 'n/a'):.2f} | Kurtosis = {row.get('kurtosis', 'n/a'):.2f} | Zero-crossing noise similarity = {zc_str}")
        # logging.info(f"    Normality test passed: {row.get('normality_pass', 'n/a')}")
        # logging.info(f"    Residuals: Ljung-Box [{ljung}], {row.get('bg_test', 'BG?')} [{bg}], Normality [{norm}]")
        # logging.info(f"    Custom Corr [{custom_corr_symbol}]")
    
        # comp_stats = row.get("composite_stats", {})
        # if comp_stats:
            # logging.info(f"    Composite stats: "
                        #  f"Pearson = {comp_stats['pearson_corr']:.2f}, "
                        #  f"Spearman = {comp_stats['spearman_corr']:.2f}, "
                        #  f"Spectral = {comp_stats['spectral_ratio']:.2f}, "
                        #  f"R² = {comp_stats['avg_rolling_r2']:.2f}, "
                        #  f"Run ratio = {comp_stats['run_ratio']:.2f}"
                        # )
        
    return sorted_models

def advanced_residual_diagnostics(H0, residuals, model_name, enable_tests=True, enable_custom_corr=True):
    if not enable_tests:
        logging.info(f"Skipping residual diagnostics for {model_name} due to CLI flag.")
        return {
            "skewness": None,
            "kurtosis": None,
            "normality_pass": None,
            "ljung_stat": None,
            "ljung_p": None,
            "ljung_failed": None,
            # "bg_test": None,
            # "bg_stat": None,
            # "bg_p": None,
            # "bg_failed": None,
            "reset_stat": None,
            "reset_p": None,
            "reset_failed": None,
            # "cooks_max": None,
            # "cooks_extreme": None,
            "zero_crossing_similarity": None
        }

    autocorr = autocorrelation_tests(H0, residuals, model_name)

    residuals_std = (residuals - np.mean(residuals)) / np.std(residuals)
    normality_pass = np.all(np.abs(residuals_std) < 3)

    skewness = pd.Series(residuals).skew()
    kurtosis = pd.Series(residuals).kurtosis()

    # logging.info(f"Additional diagnostics for {model_name}:")
    # logging.info(f"Skewness: {skewness:.3f}, Kurtosis: {kurtosis:.3f}")
    # if not normality_pass:
    #     logging.warning("Residuals may not be normally distributed (outliers or heavy tails).")

    result = {
        **autocorr,
        "skewness": skewness,
        "kurtosis": kurtosis,
        "normality_pass": normality_pass
    }

    if enable_custom_corr:
        custom = custom_residual_pattern_test(residuals)
        result.update({
            "composite_flagged": custom.get("composite_flagged"),
            "composite_stats": custom.get("composite_stats")
        })

        # stats = custom["composite_stats"]
        # if stats:
            # logging.info("Composite residual correlation test:")
            # logging.info(f"  Lag-1 Pearson = {stats['pearson_corr']:.3f}")
            # logging.info(f"  Lag-1 Spearman = {stats['spearman_corr']:.3f}")
            # logging.info(f"  Spectral ratio = {stats['spectral_ratio']:.3f}")
            # logging.info(f"  Avg rolling R² = {stats['avg_rolling_r2']:.3f}")
            # logging.info(f"  Run ratio = {stats['run_ratio']:.3f}")

        # if custom["composite_flagged"]:
        #     logging.warning("⚠️ Residuals flagged as correlated by composite test.")
        # else:
        #     logging.info("✓ Composite test: residuals appear uncorrelated.")

    return result

def process_csv_files_in_folder(config, skip_tests=False, plot_normalized=False):
    input_folder = config["general"]["input_dir"]
    output_folder = config["general"]["results_dir"]
    maxfev = config["general"]["maxfev"]
    lags = config["general"]["lags"]

    skip_tests = config.get("cli_flags", {}).get("skip_tests", False)
    skip_normres = config.get("cli_flags", {}).get("no_normalized", False)
    custom_corr_enabled = config.get("cli_flags", {}).get("custom_residual_check", True)

    global_delta_range = collect_global_max_deltadelta(input_folder)
    if global_delta_range is None:
        logging.warning("Global Δδ range could not be determined. Normalized residuals may be inaccurate.")

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
                logging.info(f"\nEvaluating model: {model_name}")
                guess = smart_initial_guess(model['lambda'], model['initial_guess'], model['bounds'], H0, d_delta_exp)
                bounds = model['bounds']
                func = model['lambda']

                params, cov = curve_fit(func, H0, d_delta_exp, p0=guess, bounds=bounds, maxfev=maxfev)
                fit_vals = func(H0, *params)
                residuals = fit_vals - d_delta_exp

                delta_max = np.max(np.abs(d_delta_exp)) if np.max(np.abs(d_delta_exp)) > 0 else 1
                
                normalized_residuals = residuals / global_delta_range
                weighted_rmse = np.sqrt(np.mean(normalized_residuals ** 2))
                
                

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

                diagnostics = advanced_residual_diagnostics(
                    H0,
                    residuals,
                    model_name,
                    enable_tests=not skip_tests,
                    enable_custom_corr=True  
                )

                output_rows.append({
                    "file": filename,
                    "model": model_name,
                    "parameters": params.tolist(),
                    "standard_errors": std_err.tolist(),
                    "r_squared": r2,
                    "AIC": aic,
                    "BIC": bic,
                    "RMSE": rmse,
                    "weighted_RMSE": weighted_rmse,
                    "confidence_intervals": [(v - 1.96*s, v + 1.96*s) for v, s in zip(params, std_err)],
                    "fitted_values": fit_vals.tolist(),
                    "residuals": residuals.tolist(),
                    "normalized_residuals": normalized_residuals if not skip_normres else [],
                    "H_over_G": (H0 / G0).tolist(),
                    **diagnostics
                })

                logging.info(f"Model {model_name} fit completed")
                logging.info(f"Parameters: {params}")
                # logging.info(f"R²: {r2:.4f}, AIC: {aic:.2f}, BIC: {bic:.2f}, RMSE: {rmse:.4f}, Weighted RMSE: {weighted_rmse:.4f}")

            except Exception as e:
                logging.error(f"Exception in model {model_name}: {e}")
                logging.error(traceback.format_exc())
                continue

        sorted_models = compare_models_by_metric(output_rows, metric="AIC")

        output_file = os.path.join(output_folder, filename.replace(".csv", "_results.csv"))
        save_combined_csv(sorted_models, output_file)

        plot_path = os.path.join(output_folder, filename.replace(".csv", "_plot.png"))
        plot_results(H0, G0, d_delta_exp, results, plot_path, file_title=filename, show_normalized=not skip_normres)

def plot_results(H0, G0, d_delta_exp, model_results, filename, file_title=None, show_normalized=True):
    fig, axes = plt.subplots(nrows=5, ncols=3 if show_normalized else 2, figsize=(18, 14))
    model_names = list(model_results.keys())

    if file_title:
        fig.suptitle(f"File: {file_title}", fontsize=18, fontweight='bold', color='black', y=1.02)

    for i, model in enumerate(model_names):
        fitted, residuals = model_results[model]
        x = H0 / G0

        delta_max = np.max(np.abs(d_delta_exp))
        norm_residuals = residuals / delta_max if delta_max > 0 else residuals

        # Fit
        axes[i, 0].scatter(x, d_delta_exp, label='Experimental')
        axes[i, 0].plot(x, fitted, color='red', label=f'{model} Fit')
        axes[i, 0].set_title(f'{model} Fit')
        axes[i, 0].set_xlabel('[H]/[G]')
        axes[i, 0].set_ylabel('Δδ [Hz]')
        axes[i, 0].legend()

        # Residuals
        axes[i, 1].scatter(x, residuals, color='red')
        axes[i, 1].axhline(0, linestyle='--')
        axes[i, 1].set_title(f'{model} Residuals')
        axes[i, 1].set_xlabel('[H]/[G]')
        axes[i, 1].set_ylabel('Residual [Hz]')

        # Normalized Residuals (if enabled)
        if show_normalized:
            axes[i, 2].scatter(x, norm_residuals, color='purple')
            axes[i, 2].axhline(0, linestyle='--')
            axes[i, 2].set_title(f'{model} Normalized Residuals')
            axes[i, 2].set_xlabel('[H]/[G]')
            axes[i, 2].set_ylabel('Residual / Max')

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.subplots_adjust(top=0.9)
    plt.savefig(filename)
    plt.close()
