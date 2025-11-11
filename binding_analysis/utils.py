"""
Utility functions for binding analysis.

This module contains helper functions for data validation, residual analysis,
statistical tests, and file operations used in the binding analysis workflow.

Author: Filip Hládek
License: MIT
"""

import logging
import os
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.fft import rfft
from scipy.optimize import curve_fit
from scipy.stats import kurtosis, normaltest, pearsonr, skew, spearmanr
from statsmodels.sandbox.stats.runs import runstest_1samp
from statsmodels.stats.diagnostic import (
    acorr_breusch_godfrey,
    acorr_ljungbox,
    het_white,
)

# Type aliases
ArrayLike = Union[np.ndarray, list, float]


def custom_residual_pattern_test(residuals: ArrayLike) -> Dict[str, Any]:
    """
    Perform custom pattern-based tests on residuals.

    Analyzes residuals for systematic patterns using multiple statistical tests
    including correlation, spectral analysis, rolling R², and runs tests.

    Args:
        residuals: Array of residual values

    Returns:
        Dictionary containing test statistics and overall flag indicating
        if systematic patterns were detected
    """
    residuals = np.asarray(residuals)

    if len(residuals) < 9:
        return {
            "composite_stats": None,
            "composite_flagged": None,
            "message": "Insufficient data points for pattern analysis",
        }

    n = len(residuals)

    # Lag-1 Pearson and Spearman correlation tests
    try:
        pearson_corr, _ = pearsonr(residuals[:-1], residuals[1:])
        spearman_corr, _ = spearmanr(residuals[:-1], residuals[1:])
    except Exception:
        pearson_corr = spearman_corr = 0.0

    # Spectral energy ratio (low frequency vs total)
    try:
        freqs = rfft(residuals)
        energy = np.abs(freqs) ** 2
        low_freq_energy = np.sum(energy[: max(1, n // 10)])
        total_energy = np.sum(energy)
        spectral_ratio = low_freq_energy / total_energy if total_energy > 0 else 0
    except Exception:
        spectral_ratio = 0.0

    # Rolling R² (linear fit on residual windows)
    window = min(5, len(residuals) - 2)
    rolling_r2 = []

    if window > 2:
        for i in range(n - window + 1):
            x = np.arange(window)
            y = residuals[i : i + window]
            try:
                slope, intercept = np.polyfit(x, y, 1)
                y_pred = slope * x + intercept
                ss_res = np.sum((y - y_pred) ** 2)
                ss_tot = np.sum((y - np.mean(y)) ** 2)
                r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
                rolling_r2.append(r2)
            except Exception:
                rolling_r2.append(0.0)

    avg_r2 = np.mean(rolling_r2) if rolling_r2 else 0.0

    # Runs test (sign changes)
    try:
        signs = np.sign(residuals - np.mean(residuals))
        run_stat, p_run = runstest_1samp(signs)
        expected_runs = (2 * np.sum(signs != 0) - 1) / 3
        actual_runs = np.sum(np.diff(signs) != 0)
        run_ratio = actual_runs / expected_runs if expected_runs > 0 else 1.0
    except Exception:
        run_ratio = 1.0

    # Heuristic thresholds for pattern detection
    flagged = (
        abs(pearson_corr) > 0.35
        or abs(spearman_corr) > 0.35
        or spectral_ratio > 0.3
        or avg_r2 > 0.35
        or run_ratio < 0.65
        or run_ratio > 1.35
    )

    return {
        "composite_stats": {
            "pearson_corr": pearson_corr,
            "spearman_corr": spearman_corr,
            "spectral_ratio": spectral_ratio,
            "avg_rolling_r2": avg_r2,
            "run_ratio": run_ratio,
        },
        "composite_flagged": flagged,
    }


def collect_global_max_deltadelta(input_folder: str) -> float:
    """
    Scan all CSV files in the input folder and return the maximum Δδ across all files.

    This function is used for global normalization of residuals when the
    --global_norm flag is used.

    Args:
        input_folder: Path to directory containing CSV files

    Returns:
        Maximum absolute chemical shift change across all files
    """
    input_path = Path(input_folder)
    max_ddeltas = []

    for csv_file in input_path.glob("*.csv"):
        try:
            df = pd.read_csv(csv_file)
            if "delta" not in df.columns:
                logging.warning(f"No 'delta' column found in {csv_file.name}")
                continue

            # Calculate delta-delta (change from first point)
            delta_changes = np.abs(df["delta"] - df["delta"].iloc[0])
            delta_changes.iloc[0] = 0  # First point is always 0
            max_ddeltas.append(delta_changes.max())

        except Exception as e:
            logging.warning(f"Failed to parse {csv_file.name}: {e}")

    if not max_ddeltas:
        logging.warning("No valid CSV files found or processed")
        return 1.0

    return max(max_ddeltas)


def smart_initial_guess(
    model_func: Callable,
    H0: ArrayLike,
    d_delta_exp: ArrayLike,
    default_guess: List[float],
    bounds: Tuple[List[float], List[float]],
    max_iter: int = 5,
    perturbation: float = 0.1,
) -> List[float]:
    """
    Generate intelligent initial parameter guesses by testing variations.

    This function tries multiple perturbed versions of the default guess
    and returns the one that gives the best initial fit.

    Args:
        model_func: Model function to test
        H0: Host concentration data
        d_delta_exp: Experimental data to fit
        default_guess: Default initial parameter guess
        bounds: Parameter bounds (lower, upper)
        max_iter: Maximum number of guessing attempts
        perturbation: Relative perturbation size (fraction)

    Returns:
        Best initial guess found
    """
    best_guess = list(default_guess)
    best_score = np.inf

    for _ in range(max_iter):
        # Generate perturbed guess
        guess = [
            p * (1 + np.random.uniform(-perturbation, perturbation))
            for p in default_guess
        ]

        # Ensure guess respects bounds
        lower_bounds, upper_bounds = bounds
        for i in range(len(guess)):
            guess[i] = np.clip(guess[i], lower_bounds[i], upper_bounds[i])

        try:
            popt, _ = curve_fit(
                model_func, H0, d_delta_exp, p0=guess, bounds=bounds, maxfev=100000
            )
            residuals = model_func(H0, *popt) - d_delta_exp
            score = np.sum(residuals**2)

            if score < best_score:
                best_score = score
                best_guess = list(guess)

        except Exception:
            continue

    return best_guess


def compare_models_by_metric(
    output_rows: List[Dict[str, Any]], metric: str = "AIC"
) -> None:
    """
    Compare and rank models based on a specified metric.

    Args:
        output_rows: List of model result dictionaries
        metric: Metric to use for comparison (default: "AIC")
    """
    try:
        model_metrics = [
            (row["model"], row[metric]) for row in output_rows if metric in row
        ]
        if not model_metrics:
            logging.warning(f"No models contain the metric '{metric}'")
            return

        sorted_models = sorted(model_metrics, key=lambda x: x[1])
        logging.info(f"Model Ranking by {metric} (lower is better):")
        for rank, (model, val) in enumerate(sorted_models, start=1):
            logging.info(f"{rank}. {model}: {metric} = {val:.2f}")
    except Exception as e:
        logging.warning(f"Model comparison skipped due to: {e}")


def advanced_residual_diagnostics(
    residuals: ArrayLike, model_name: str
) -> Dict[str, Any]:
    """
    Perform advanced statistical analysis on residuals.

    Calculates skewness, kurtosis, and normality tests to assess
    the quality of the residual distribution.

    Args:
        residuals: Array of residual values
        model_name: Name of the model for logging

    Returns:
        Dictionary containing diagnostic statistics
    """
    try:
        residuals = np.asarray(residuals)
        s = skew(residuals)
        k = kurtosis(residuals)
        stat, p = normaltest(residuals)

        return {"skew": s, "kurtosis": k, "normality_p": p, "normality_pass": p > 0.05}
    except Exception as e:
        logging.warning(f"Residual diagnostics for {model_name} skipped due to: {e}")
        return {}


def validate_data(df: pd.DataFrame) -> bool:
    """
    Validate that the input dataframe contains required columns and data.

    Args:
        df: Input dataframe to validate

    Returns:
        True if data is valid, False otherwise
    """
    required_cols = ["H", "G", "delta"]
    missing = [col for col in required_cols if col not in df.columns]

    if missing:
        logging.error(f"Missing required columns: {missing}")
        return False

    if df[required_cols].isnull().any().any():
        logging.warning("Detected NaN values in input data")

    if len(df) < 3:
        logging.error("Insufficient data points (minimum 3 required)")
        return False

    return True


def autocorrelation_tests(H0, residuals, model_name, lags=10):
    results = {}
    n = len(H0)
    X = sm.add_constant(H0)

    lags = min(10, n // 5)

    # Ljung-Box
    if n >= 2 * lags + 1:
        lb_test = acorr_ljungbox(residuals, lags=[min(lags, n - 2)], return_df=True)
        lb_stat = lb_test["lb_stat"].values[0]
        lb_p = lb_test["lb_pvalue"].values[0]
        results.update(
            {"ljung_stat": lb_stat, "ljung_p": lb_p, "ljung_failed": lb_p < 0.05}
        )
    else:
        logging.info("Ljung-Box test skipped (too few data points).")
        results.update({"ljung_stat": None, "ljung_p": None, "ljung_failed": None})

    # Breusch-Godfrey or White
    try:
        model = sm.OLS(residuals, X).fit()
        bg_stat, bg_p, _, _ = acorr_breusch_godfrey(model, nlags=min(lags, n - 2))
        bg_name = "Breusch-Godfrey"
    except Exception:
        model = sm.OLS(residuals, X).fit()
        bg_stat, bg_p, _, _ = het_white(model.resid, model.model.exog)
        bg_name = "White's Test"

    results.update(
        {"bg_test": bg_name, "bg_stat": bg_stat, "bg_p": bg_p, "bg_failed": bg_p < 0.05}
    )

    # Ramsey RESET test
    if n >= 15:
        try:
            model = sm.OLS(residuals, X).fit()
            from statsmodels.stats.diagnostic import linear_reset

            reset_test = linear_reset(model, power=2, use_f=True)
            results.update(
                {
                    "reset_stat": reset_test.fvalue,
                    "reset_p": reset_test.pvalue,
                    "reset_failed": reset_test.pvalue < 0.05,
                }
            )
        except Exception as e:
            logging.warning(f"RESET test failed: {e}")
    else:
        logging.info("Ramsey RESET test skipped (too few data points).")

    # Cook's distance
    if n >= 10:
        try:
            influence = model.get_influence()
            cooks_d = influence.cooks_distance[0]
            max_cook = np.max(cooks_d)
            n_extreme = np.sum(cooks_d > (4 / n))
            results.update({"cooks_max": max_cook, "cooks_extreme": n_extreme})
            # Add logging to make Cook's Distance visible in output
            logging.info(f"Cook's Distance: max = {max_cook:.4f}, extreme (>4/n): {n_extreme}")
        except Exception as e:
            logging.warning(f"Cook's Distance failed: {e}")
    else:
        logging.info("Cook's Distance skipped (too few data points).")

    # Zero-crossing randomness test
    try:
        zero_crossings = np.sum(np.diff(np.sign(residuals)) != 0)
        simulated_crossings = []
        for _ in range(1000):
            sim = np.random.normal(0, 1, n)
            simulated_crossings.append(np.sum(np.diff(np.sign(sim)) != 0))
        sim_mean = np.mean(simulated_crossings)
        similarity = 100 * (1 - abs(zero_crossings - sim_mean) / sim_mean)
        results["zero_crossings"] = zero_crossings
        results["crossing_similarity"] = similarity
        # Add logging to make zero-crossings visible in output
        logging.info(f"Zero-crossings: {zero_crossings} | Similarity to white noise: {similarity:.2f}%")
    except Exception as e:
        logging.warning(f"Zero-crossing similarity test failed: {e}")

    return results


def delete_old_result_files(folder_path: str) -> None:
    """
    Clean up old result files from previous runs.

    Args:
        folder_path: Path to directory to clean
    """
    folder = Path(folder_path)
    if not folder.exists():
        return

    for f in folder.iterdir():
        if f.suffix in [".csv", ".png"] and ("_results" in f.name or "_plot" in f.name):
            try:
                f.unlink()
                logging.info(f"Deleted: {f.name}")
            except Exception as e:
                logging.warning(f"Failed to delete {f.name}: {e}")


def save_combined_csv(results: List[Dict[str, Any]], output_file: str) -> None:
    """
    Save analysis results to a combined CSV file.

    This function consolidates all model fitting results, residuals, and
    diagnostic statistics into a comprehensive CSV file for further analysis.

    Args:
        results: List of result dictionaries from model fitting
        output_file: Path to output CSV file
    """
    rows = []

    for result in results:
        num_points = len(result["fitted_values"])

        # Add data points and residuals
        for i in range(num_points):
            row = {
                "file": result["file"],
                "model": result["model"],
                "H_over_G": result["H_over_G"][i],
                "fitted_value": result["fitted_values"][i],
                "residual": result["residuals"][i],
            }
            if result.get("normalized_residuals"):
                row["normalized_residual"] = result["normalized_residuals"][i]
            rows.append(row)

        # Add parameter information
        param_names = result.get(
            "parameter_names",
            [f"param_{i+1}" for i in range(len(result["parameters"]))],
        )
        for i, param in enumerate(result["parameters"]):
            param_name = param_names[i] if i < len(param_names) else f"param_{i+1}"
            rows.append(
                {
                    "file": result["file"],
                    "model": result["model"],
                    "parameter_name": param_name,
                    "parameter_index": i + 1,
                    "parameter_value": param,
                    "standard_error": (
                        result["standard_errors"][i]
                        if i < len(result["standard_errors"])
                        else None
                    ),
                    "confidence_interval_lower": (
                        result["confidence_intervals"][i][0]
                        if i < len(result["confidence_intervals"])
                        else None
                    ),
                    "confidence_interval_upper": (
                        result["confidence_intervals"][i][1]
                        if i < len(result["confidence_intervals"])
                        else None
                    ),
                    "r_squared": result["r_squared"],
                    "AIC": result["AIC"],
                    "BIC": result["BIC"],
                    "RMSE": result["RMSE"],
                    "weighted_RMSE": result.get("weighted_RMSE"),
                    "ljung_stat": result.get("ljung_stat"),
                    "ljung_p": result.get("ljung_p"),
                    "ljung_failed": result.get("ljung_failed"),
                    "bg_test": result.get("bg_test"),
                    "bg_stat": result.get("bg_stat"),
                    "bg_p": result.get("bg_p"),
                    "bg_failed": result.get("bg_failed"),
                    "reset_stat": result.get("reset_stat"),
                    "reset_p": result.get("reset_p"),
                    "reset_failed": result.get("reset_failed"),
                    "cooks_max": result.get("cooks_max"),
                    "cooks_extreme": result.get("cooks_extreme"),
                    "crossing_similarity": result.get("crossing_similarity"),
                    "zero_crossings": result.get("zero_crossings"),
                    "composite_flagged": result.get("composite_flagged"),
                    "composite_stats": str(result.get("composite_stats")),
                    # Additional diagnostic fields that were missing
                    "normality_pass": result.get("normality_pass"),
                    "shapiro_stat": result.get("shapiro_stat"),
                    "shapiro_p": result.get("shapiro_p"),
                    "skewness": result.get("skewness"),
                    "kurtosis": result.get("kurtosis"),
                    "pearson_corr": result.get("composite_stats", {}).get("pearson_corr") if result.get("composite_stats") else None,
                    "spearman_corr": result.get("composite_stats", {}).get("spearman_corr") if result.get("composite_stats") else None, 
                    "spectral_ratio": result.get("composite_stats", {}).get("spectral_ratio") if result.get("composite_stats") else None,
                    "avg_rolling_r2": result.get("composite_stats", {}).get("avg_rolling_r2") if result.get("composite_stats") else None,
                    "run_ratio": result.get("composite_stats", {}).get("run_ratio") if result.get("composite_stats") else None,
                }
            )

    df = pd.DataFrame(rows)
    df.to_csv(output_file, index=False)
    logging.info(f"Combined results saved to: {output_file}")
