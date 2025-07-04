"""
Main analysis module for binding isotherm fitting.

This module contains the core analysis functions for fitting NMR titration data
to various thermodynamic binding models and performing statistical diagnostics.

Author: Filip Hládek
License: MIT
"""

import logging
import os
import traceback
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from .models import model_definitions
from scipy.optimize import curve_fit, least_squares
from scipy.stats import pearsonr, spearmanr
from statsmodels.api import OLS, add_constant
from tabulate import tabulate
from .utils import (
    autocorrelation_tests,
    collect_global_max_deltadelta,
    custom_residual_pattern_test,
    save_combined_csv,
    validate_data,
)

# Type aliases for clarity
ArrayLike = Union[np.ndarray, list, float]


def smart_initial_guess(
    model_func,
    guess: List[float],
    bounds: Tuple[List[float], List[float]],
    H0: np.ndarray,
    d_delta_exp: np.ndarray,
    step: int = 10,
    max_iter: int = 10,
) -> np.ndarray:
    """
    Intelligently refine initial parameter guesses through grid search.

    This function systematically perturbs each parameter in the initial guess
    to find a better starting point for the curve fitting.

    Args:
        model_func: The model function to test
        guess: Initial parameter guess
        bounds: Parameter bounds (lower, upper)
        H0: Host concentration data
        d_delta_exp: Experimental data to fit
        step: Step size for parameter perturbation
        max_iter: Maximum number of iterations

    Returns:
        Improved initial guess
    """
    best_guess = np.array(guess)
    best_rmse = np.inf

    for iteration in range(max_iter):
        improved = False

        for i in range(len(best_guess)):
            for delta in [-step, step]:
                candidate = best_guess.copy()
                candidate[i] += delta

                # Ensure bounds are respected
                lower_bounds, upper_bounds = bounds
                candidate[i] = np.clip(candidate[i], lower_bounds[i], upper_bounds[i])

                try:
                    fit_vals = model_func(H0, *candidate)
                    rmse = np.sqrt(np.mean((fit_vals - d_delta_exp) ** 2))

                    if rmse < best_rmse:
                        best_guess = candidate
                        best_rmse = rmse
                        improved = True

                except Exception:
                    # Skip if model evaluation fails
                    continue

        # Early exit if no improvement found
        if not improved:
            break

    return best_guess


def metric_value_range(metric: str) -> str:
    """
    Get the expected value range for a given metric.

    Args:
        metric: Name of the metric

    Returns:
        String describing the expected range
    """
    ranges = {
        "r2": "0 to 1",
        "rmse": "0 to ∞",
        "wrmse": "0 to ∞",
        "aic": "−∞ to ∞",
        "bic": "−∞ to ∞",
        "normality": "✓ / ✗",
        "ljung": "0 to 1 (p-value)",
        "reset": "0 to 1 (p-value)",
        "pearson": "−1 to 1",
        "spearman": "−1 to 1",
        "rolling_r2": "0 to 1",
    }
    return ranges.get(metric, "n/a")


def interpret_diagnostic(
    metric: str, value: float, threshold: float, passed: bool
) -> str:
    """
    Provide interpretation for diagnostic test results.

    Args:
        metric: Name of the diagnostic metric
        value: Test statistic value
        threshold: Significance threshold
        passed: Whether the test passed

    Returns:
        Interpretation string explaining the result
    """
    if passed:
        return "OK - No significant issues detected"

    explanations = {
        "normality": (
            "Residuals are not normally distributed. This could mean outliers, skewed error structure, "
            "or heavier tails than expected. It may affect confidence in fitted parameters."
        ),
        "ljung": (
            "Residuals show repeating patterns or cycles, "
            "indicating poor model fit over time or concentration range."
        ),
        "reset": (
            "Model structure may be too simple. There may be a missing non-linear term or interaction "
            "not captured by the current model equation."
        ),
        "pearson": (
            "Residuals increase or decrease in a linear fashion across the domain. Suggests the model "
            "is systematically over- or under-predicting."
        ),
        "spearman": (
            "Residuals follow a consistent monotonic trend (e.g., always increasing or decreasing). "
            "Implies model misses a non-linear or saturating effect."
        ),
        "rolling_r2": (
            "Local trends or structure detected. Residuals show 'blocks' or gradual shifts, "
            "possibly due to unmodeled system behavior or missing covariates."
        ),
        "spectral_ratio": (
            "Residuals contain dominant low-frequency components, possibly indicating periodic behavior "
            "or systematic deviations."
        ),
        "run_ratio": (
            "Number of sign changes in residuals deviates from expected. Indicates clustering or directional drift."
        ),
        "zero_crossings": (
            "Residuals show fewer random fluctuations than expected. White-noise similarity is low, "
            "suggesting model misfit or trend."
        ),
        "cooks_distance": (
            "One or more data points may exert excessive influence on the model fit. Indicates potential outliers."
        ),
        "bg_p": (
            "Higher-order autocorrelation detected. Residuals are not independent across lags, "
            "suggesting temporal or concentration structure not captured by the model."
        ),
    }

    return explanations.get(
        metric, "Check the residuals for patterns not explained by the model."
    )


def compare_models_by_metric(
    output_rows: List[Dict[str, Any]], metric: str = "AIC"
) -> None:
    """
    Compare and rank models based on a specified metric.

    Args:
        output_rows: List of model result dictionaries
        metric: Metric to use for comparison (default: "AIC")
    """
    if not output_rows:
        logging.warning("No models to compare")
        return

    # Filter rows that have the requested metric
    valid_rows = [
        row for row in output_rows if metric in row and row[metric] is not None
    ]

    if not valid_rows:
        logging.warning(f"No models contain valid '{metric}' values")
        return

    sorted_models = sorted(valid_rows, key=lambda r: r[metric])
    logging.info(f"\nModel ranking by {metric} (lower is better):")

    for rank, row in enumerate(sorted_models, 1):
        model = row["model"]
        aic = row.get("AIC", "N/A")
        bic = row.get("BIC", "N/A")
        r2 = row.get("r_squared", "N/A")
        rmse = row.get("RMSE", "N/A")
        wrmse = row.get("weighted_RMSE", "N/A")

        # Format metric value
        metric_val = row[metric]
        if isinstance(metric_val, (int, float)):
            metric_str = f"{metric_val:.3f}"
        else:
            metric_str = str(metric_val)

        logging.info(f"{rank}. {model}: {metric} = {metric_str}")
        logging.info(
            f"    R² = {r2:.4f}, RMSE = {rmse:.4f}"
            if isinstance(r2, (int, float)) and isinstance(rmse, (int, float))
            else ""
        )

        # Report any diagnostic flags
        flags = []
        if row.get("ljung_failed"):
            flags.append("Autocorr")
        if row.get("reset_failed"):
            flags.append("Reset")
        if row.get("composite_flagged"):
            flags.append("Pattern")

        if flags:
            logging.info(f"    Diagnostic flags: {', '.join(flags)}")

        logging.info("-" * 50)
        # Extract diagnostic values from row data
        cooks = row.get("cooks_max")
        zero_crossings = row.get("zero_crossings")
        crossing_sim = row.get("crossing_similarity")
        bg_white = row.get("bg_p")
        bg_or_white = row.get("bg_test")
        ljung = row.get("ljung_p")
        reset_p = row.get("reset_p")
        spectral = row.get("spectral_ratio")
        run_ratio = row.get("run_ratio")
        comp_stats = row.get("composite_stats", {})

        # Build section headers and metric keys
        section_headers = {
            "Criteria": ["R²", "RMSE", "Weighted RMSE", "AIC", "BIC"],
            "Main Tests": [
                "Normality test",
                "Ljung-Box p",
                "RESET p",
                "Pearson corr",
                "Spearman corr",
                "Rolling R²",
            ],
            "Optional Tests": [
                "Spectral ratio",
                "Run ratio",
                "Zero-crossings",
                "Cook’s Distance",
                "Breusch-Godfrey",
            ],
        }

        table_data = [
            ["Metric", "Value", "Acceptable Range", "Possible Range", "Interpretation"]
        ]

        for section, metrics in section_headers.items():
            full_span = f"=== {section} ===".center(95)
            table_data.append([full_span, "", "", "", ""])

            for metric in metrics:
                if metric == "R²":
                    table_data.append(
                        [
                            metric,
                            f"{r2:.4f}",
                            "> 0.98",
                            "0 to 1",
                            interpret_diagnostic("r2", r2, 0.98, r2 > 0.98),
                        ]
                    )
                elif metric == "RMSE":
                    table_data.append(
                        [
                            metric,
                            f"{rmse:.4f}",
                            "Low",
                            "0 to ∞",
                            "Lower indicates better fit",
                        ]
                    )
                elif metric == "Weighted RMSE" and wrmse is not None:
                    table_data.append(
                        [
                            metric,
                            f"{wrmse:.4f}",
                            "Low",
                            "0 to ∞",
                            "Used for comparison across datasets",
                        ]
                    )
                elif metric == "AIC":
                    table_data.append(
                        [
                            metric,
                            f"{aic:.2f}",
                            "Low",
                            "−∞ to ∞",
                            "Used to compare model fit (penalized)",
                        ]
                    )
                elif metric == "BIC":
                    table_data.append(
                        [
                            metric,
                            f"{bic:.2f}",
                            "Low",
                            "−∞ to ∞",
                            "Stricter than AIC for model comparison",
                        ]
                    )
                elif metric == "Normality test":
                    normality_pass = row.get("normality_pass", True)
                    table_data.append(
                        [
                            metric,
                            "Passed" if normality_pass else "Failed",
                            "Passed",
                            "Passed / Failed",
                            interpret_diagnostic("normality", None, None, normality_pass),
                        ]
                    )
                elif metric == "Ljung-Box p" and ljung is not None:
                    table_data.append(
                        [
                            metric,
                            f"{ljung:.3f}",
                            "> 0.05",
                            "0 to 1",
                            interpret_diagnostic("ljung", ljung, 0.05, ljung > 0.05),
                        ]
                    )
                elif metric == "RESET p" and reset_p is not None:
                    table_data.append(
                        [
                            metric,
                            f"{reset_p:.4f}",
                            "> 0.05",
                            "0 to 1",
                            interpret_diagnostic(
                                "reset", reset_p, 0.05, reset_p > 0.05
                            ),
                        ]
                    )
                elif metric == "Pearson corr" and "pearson_corr" in comp_stats:
                    p = comp_stats["pearson_corr"]
                    table_data.append(
                        [
                            metric,
                            f"{p:.2f}",
                            "< 0.35",
                            "-1 to 1",
                            interpret_diagnostic("pearson", p, 0.35, abs(p) < 0.35),
                        ]
                    )
                elif metric == "Spearman corr" and "spearman_corr" in comp_stats:
                    s = comp_stats["spearman_corr"]
                    table_data.append(
                        [
                            metric,
                            f"{s:.2f}",
                            "< 0.35",
                            "-1 to 1",
                            interpret_diagnostic("spearman", s, 0.35, abs(s) < 0.35),
                        ]
                    )
                elif metric == "Rolling R²" and "avg_rolling_r2" in comp_stats:
                    r = comp_stats["avg_rolling_r2"]
                    table_data.append(
                        [
                            metric,
                            f"{r:.2f}",
                            "< 0.35",
                            "0 to 1",
                            interpret_diagnostic("rolling_r2", r, 0.35, r < 0.35),
                        ]
                    )
                elif metric == "Spectral ratio" and spectral is not None:
                    table_data.append(
                        [
                            metric,
                            f"{spectral:.2f}",
                            "< 0.3",
                            "0 to 1",
                            interpret_diagnostic(
                                "spectral_ratio", spectral, 0.3, spectral < 0.3
                            ),
                        ]
                    )
                elif metric == "Run ratio" and run_ratio is not None:
                    passed = 0.65 < run_ratio < 1.35
                    table_data.append(
                        [
                            metric,
                            f"{run_ratio:.2f}",
                            "≈ 1.0",
                            "0 to ~2.0",
                            interpret_diagnostic("run_ratio", run_ratio, None, passed),
                        ]
                    )
                elif metric == "Zero-crossings" and crossing_sim is not None:
                    passed = crossing_sim > 80
                    table_data.append(
                        [
                            metric,
                            f"{crossing_sim:.1f}%",
                            "> 80%",
                            "0 to 100%",
                            interpret_diagnostic(
                                "zero_crossings", crossing_sim, 80, passed
                            ),
                        ]
                    )

                elif metric == "Cook’s Distance" and cooks is not None:
                    passed = cooks < 1.0
                    table_data.append(
                        [
                            metric,
                            f"{cooks:.3f}",
                            "< 1.0",
                            "0 to ∞",
                            interpret_diagnostic("cooks_distance", cooks, 1.0, passed),
                        ]
                    )

                elif metric == "Breusch-Godfrey" and bg_white is not None:
                    passed = bg_white > 0.05
                    table_data.append(
                        [
                            bg_or_white,
                            f"{bg_white:.3f}",
                            "> 0.05",
                            "0 to 1",
                            interpret_diagnostic("bg_p", bg_white, 0.05, passed),
                        ]
                    )

        logging.info(
            f"\n{rank}. Model: {model}\n"
            + tabulate(table_data, headers="firstrow", tablefmt="fancy_grid")
        )

    return sorted_models


def advanced_residual_diagnostics(
    H0: ArrayLike,
    residuals: ArrayLike,
    model_name: str,
    enable_tests: bool = True,
    enable_custom_corr: bool = True,
) -> Dict[str, Any]:
    """
    Perform comprehensive residual diagnostics.

    Args:
        H0: Host concentration data
        residuals: Model residuals
        model_name: Name of the model
        enable_tests: Whether to run statistical tests
        enable_custom_corr: Whether to run custom correlation tests

    Returns:
        Dictionary containing all diagnostic results
    """
    if not enable_tests:
        logging.info(f"Skipping residual diagnostics for {model_name} due to CLI flag.")
        return {
            "skewness": None,
            "kurtosis": None,
            "normality_pass": None,
            "ljung_stat": None,
            "ljung_p": None,
            "ljung_failed": None,
            "bg_test": None,
            "bg_stat": None,
            "bg_p": None,
            "bg_failed": None,
            "reset_stat": None,
            "reset_p": None,
            "reset_failed": None,
            "cooks_max": None,
            "cooks_extreme": None,
            "zero_crossing_similarity": None,
        }

    # Run autocorrelation tests
    autocorr = autocorrelation_tests(H0, residuals, model_name)

    # Basic normality check (3-sigma rule)
    residuals_std = (residuals - np.mean(residuals)) / np.std(residuals)
    normality_pass = np.all(np.abs(residuals_std) < 3)

    # Calculate skewness and kurtosis
    import pandas as pd

    skewness = pd.Series(residuals).skew()
    kurt = pd.Series(residuals).kurtosis()

    # Run custom pattern tests if enabled
    custom_results = {}
    if enable_custom_corr:
        custom_results = custom_residual_pattern_test(residuals)

    # Combine all results
    results = {
        "skewness": skewness,
        "kurtosis": kurt,
        "normality_pass": normality_pass,
        **autocorr,
        **custom_results,
    }

    return results
    kurtosis = pd.Series(residuals).kurtosis()

    # logging.info(f"Additional diagnostics for {model_name}:")
    # logging.info(f"Skewness: {skewness:.3f}, Kurtosis: {kurtosis:.3f}")
    # if not normality_pass:
    #     logging.warning("Residuals may not be normally distributed (outliers or heavy tails).")

    result = {
        **autocorr,
        "skewness": skewness,
        "kurtosis": kurtosis,
        "normality_pass": normality_pass,
    }

    if enable_custom_corr:
        custom = custom_residual_pattern_test(residuals)
        result.update(
            {
                "composite_flagged": custom.get("composite_flagged"),
                "composite_stats": custom.get("composite_stats"),
            }
        )

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
        logging.warning(
            "Global Δδ range could not be determined. Normalized residuals may be inaccurate."
        )

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

        H0, G0, d_delta_exp = np.genfromtxt(
            full_path, delimiter=",", names=True, unpack=True
        )
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
                guess = smart_initial_guess(
                    model["lambda"],
                    model["initial_guess"],
                    model["bounds"],
                    H0,
                    d_delta_exp,
                )
                bounds = model["bounds"]
                func = model["lambda"]

                def residuals_func(params):
                    return func(H0, *params) - d_delta_exp

                result = least_squares(
                    residuals_func, x0=guess, bounds=bounds, max_nfev=maxfev
                )
                params = result.x
                fit_vals = func(H0, *params)
                residuals = fit_vals - d_delta_exp
                n_iter = result.nfev * len(H0)

                logging.info(f"Fitting completed in {n_iter} function evaluations.")

                delta_max = (
                    np.max(np.abs(d_delta_exp))
                    if np.max(np.abs(d_delta_exp)) > 0
                    else 1
                )
                normalized_residuals = residuals / global_delta_range
                weighted_rmse = np.sqrt(np.mean(normalized_residuals**2))

                # Estimating covariance (approximation)
                from scipy.linalg import svd

                _, s, VT = svd(result.jac, full_matrices=False)
                threshold = np.finfo(float).eps * max(result.jac.shape) * s[0]
                s = s[s > threshold]
                VT = VT[: s.size]
                cov = np.dot(VT.T / s**2, VT)
                std_err = np.sqrt(np.diag(cov))

                rss = np.sum(residuals**2)
                r2 = 1 - rss / np.sum((d_delta_exp - np.mean(d_delta_exp)) ** 2)
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
                    enable_custom_corr=True,
                )

                output_rows.append(
                    {
                        "file": filename,
                        "model": model_name,
                        "parameters": params.tolist(),
                        "standard_errors": std_err.tolist(),
                        "r_squared": r2,
                        "AIC": aic,
                        "BIC": bic,
                        "RMSE": rmse,
                        "weighted_RMSE": weighted_rmse,
                        "confidence_intervals": [
                            (v - 1.96 * s, v + 1.96 * s)
                            for v, s in zip(params, std_err)
                        ],
                        "fitted_values": fit_vals.tolist(),
                        "residuals": residuals.tolist(),
                        "normalized_residuals": (
                            normalized_residuals if not skip_normres else []
                        ),
                        "H_over_G": (H0 / G0).tolist(),
                        "nfev": n_iter,
                        **diagnostics,
                    }
                )

                logging.info(f"Model {model_name} fit completed")
                logging.info(f"Parameters: {params}")

            except Exception as e:
                logging.error(f"Exception in model {model_name}: {e}")
                logging.error(traceback.format_exc())
                continue

        sorted_models = compare_models_by_metric(output_rows, metric="AIC")

        output_file = os.path.join(
            output_folder, filename.replace(".csv", "_results.csv")
        )
        save_combined_csv(sorted_models, output_file)

        plot_path = os.path.join(output_folder, filename.replace(".csv", "_plot.png"))
        plot_results(
            H0,
            G0,
            d_delta_exp,
            results,
            plot_path,
            file_title=filename,
            show_normalized=not skip_normres,
        )


def plot_results(
    H0, G0, d_delta_exp, model_results, filename, file_title=None, show_normalized=True
):
    fig, axes = plt.subplots(
        nrows=5, ncols=3 if show_normalized else 2, figsize=(18, 14)
    )
    model_names = list(model_results.keys())

    if file_title:
        fig.suptitle(
            f"File: {file_title}", fontsize=18, fontweight="bold", color="black", y=1.02
        )

    for i, model in enumerate(model_names):
        fitted, residuals = model_results[model]
        x = H0 / G0

        delta_max = np.max(np.abs(d_delta_exp))
        norm_residuals = residuals / delta_max if delta_max > 0 else residuals

        # Fit
        axes[i, 0].scatter(x, d_delta_exp, label="Experimental")
        axes[i, 0].plot(x, fitted, color="red", label=f"{model} Fit")
        axes[i, 0].set_title(f"{model} Fit")
        axes[i, 0].set_xlabel("[H]/[G]")
        axes[i, 0].set_ylabel("Δδ [Hz]")
        axes[i, 0].legend()

        # Residuals
        axes[i, 1].scatter(x, residuals, color="red")
        axes[i, 1].axhline(0, linestyle="--")
        axes[i, 1].set_title(f"{model} Residuals")
        axes[i, 1].set_xlabel("[H]/[G]")
        axes[i, 1].set_ylabel("Residual [Hz]")

        # Normalized Residuals (if enabled)
        if show_normalized:
            axes[i, 2].scatter(x, norm_residuals, color="purple")
            axes[i, 2].axhline(0, linestyle="--")
            axes[i, 2].set_title(f"{model} Normalized Residuals")
            axes[i, 2].set_xlabel("[H]/[G]")
            axes[i, 2].set_ylabel("Residual / Max")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.subplots_adjust(top=0.9)
    plt.savefig(filename)
    plt.close()
