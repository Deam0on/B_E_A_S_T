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
from models import model_definitions
from scipy.optimize import curve_fit, least_squares, minimize_scalar
from scipy.stats import pearsonr, spearmanr
from statsmodels.api import OLS, add_constant
from tabulate import tabulate
from utils import (
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
    Intelligently refine initial parameter guesses using multiple strategies.

    This redesigned function combines data-driven estimation with multi-scale
    optimization to find much better starting points for curve fitting.

    Args:
        model_func: The model function to test
        guess: Initial parameter guess
        bounds: Parameter bounds (lower, upper)
        H0: Host concentration data
        d_delta_exp: Experimental data to fit
        step: Step size for parameter perturbation (legacy parameter)
        max_iter: Maximum number of iterations

    Returns:
        Improved initial guess
    """
    
    # Ensure initial guess respects bounds and has minimum values
    guess = np.array(guess, dtype=float)
    lower_bounds, upper_bounds = bounds
    
    # Apply strict bounds checking with minimum values
    for i in range(len(guess)):
        if i < len(lower_bounds) and i < len(upper_bounds):
            # Ensure Ka (usually first parameter) has minimum value > 0
            if i == 0:  # Binding constant
                min_ka = max(lower_bounds[i], 1e-6)  # Minimum Ka = 1e-6 M^-1
                guess[i] = np.clip(guess[i], min_ka, upper_bounds[i])
            else:
                guess[i] = np.clip(guess[i], lower_bounds[i], upper_bounds[i])
    
    # Phase 1: Data-driven parameter estimation
    data_driven_guess = _estimate_parameters_from_data(H0, d_delta_exp, guess, bounds)
    
    # Phase 2: Multi-scale grid search with adaptive steps
    grid_optimized_guess = _adaptive_grid_search(
        model_func, data_driven_guess, bounds, H0, d_delta_exp, max_iter
    )
    
    # Phase 3: Fine-tune with coordinate descent
    final_guess = _coordinate_descent_refinement(
        model_func, grid_optimized_guess, bounds, H0, d_delta_exp
    )
    
    # Final bounds check to ensure no invalid values
    for i in range(len(final_guess)):
        if i < len(lower_bounds) and i < len(upper_bounds):
            if i == 0:  # Binding constant
                min_ka = max(lower_bounds[i], 1e-6)
                final_guess[i] = np.clip(final_guess[i], min_ka, upper_bounds[i])
            else:
                final_guess[i] = np.clip(final_guess[i], lower_bounds[i], upper_bounds[i])
    
    return final_guess


def _estimate_parameters_from_data(H0: np.ndarray, d_delta_exp: np.ndarray, 
                                   initial_guess: List[float], 
                                   bounds: Tuple[List[float], List[float]]) -> np.ndarray:
    """
    Estimate parameters directly from experimental data characteristics.
    
    This function analyzes the binding curve shape to provide intelligent
    initial estimates for binding constants and chemical shift parameters.
    """
    guess = np.array(initial_guess, dtype=float)
    lower_bounds, upper_bounds = bounds
    
    # Basic data analysis
    max_shift = np.max(np.abs(d_delta_exp))
    min_shift = np.min(d_delta_exp)
    
    # Skip if no significant binding observed
    if max_shift < 1.0:
        return guess
    
    # Find saturation behavior (use last 20% of data points)
    n_points = len(d_delta_exp)
    tail_start = max(1, int(0.8 * n_points))
    saturation_level = np.mean(np.abs(d_delta_exp[tail_start:]))
    
    # Estimate Ka from half-saturation point
    half_sat_target = saturation_level / 2.0
    
    # Find the point closest to half-saturation
    if max_shift > half_sat_target:
        half_idx = np.argmin(np.abs(np.abs(d_delta_exp) - half_sat_target))
        half_sat_H = H0[half_idx]
        
        # Estimate Ka based on binding model assumptions
        # For 1:1 binding: Ka ≈ 1/([H]_half) when guest is in excess
        if half_sat_H > 0:
            ka_estimate = 1.0 / half_sat_H
            
            # Apply reasonable bounds based on typical NMR binding constants
            ka_estimate = np.clip(ka_estimate, 1, 10000)
            
            # Update Ka parameter (usually first parameter) with minimum bound
            if len(guess) > 0 and 0 < len(lower_bounds):
                min_ka = max(lower_bounds[0], 1e-6)  # Ensure minimum Ka > 0
                ka_bounded = np.clip(ka_estimate, min_ka, upper_bounds[0])
                guess[0] = ka_bounded
    
    # Ensure Ka is never zero or negative
    if len(guess) > 0:
        min_ka = max(lower_bounds[0] if len(lower_bounds) > 0 else 1e-6, 1e-6)
        if guess[0] <= 0:
            guess[0] = min_ka
    
    # Estimate chemical shift parameters from data magnitude
    shift_estimate = max_shift * 1.2  # Slightly larger than observed maximum
    
    # Update chemical shift parameters (usually the last parameters)
    for i in range(len(guess)):
        # Look for parameters that could be chemical shifts (unbounded below, large upper bound)
        if (i < len(lower_bounds) and 
            lower_bounds[i] == -np.inf and 
            upper_bounds[i] == np.inf):
            
            # This looks like a chemical shift parameter
            if shift_estimate != 0:
                guess[i] = shift_estimate if d_delta_exp[-1] > d_delta_exp[0] else -shift_estimate
    
    # For multi-parameter models, apply heuristics
    if len(guess) >= 4:  # Models like HG₂, H₂G
        # Second binding constant is typically weaker (statistical factor of ~4)
        if len(guess) > 1 and guess[0] > 0:
            guess[1] = max(guess[0] / 4.0, 1e-6)  # Ensure minimum value
            guess[1] = np.clip(guess[1], 
                             max(lower_bounds[1], 1e-6) if len(lower_bounds) > 1 else 1e-6, 
                             upper_bounds[1] if len(upper_bounds) > 1 else guess[0])
    
    return guess


def _adaptive_grid_search(model_func, initial_guess: np.ndarray, 
                         bounds: Tuple[List[float], List[float]],
                         H0: np.ndarray, d_delta_exp: np.ndarray,
                         max_iter: int = 10) -> np.ndarray:
    """
    Perform adaptive multi-scale grid search optimization.
    
    This uses different step sizes for different parameter types and
    progressively refines the search around the best candidates.
    """
    best_guess = initial_guess.copy()
    lower_bounds, upper_bounds = bounds
    
    # Ensure Ka parameters have minimum values
    for i in range(len(best_guess)):
        if i < len(lower_bounds):
            if i == 0 or i == 1:  # Ka and Kd parameters
                min_val = max(lower_bounds[i], 1e-6)
                best_guess[i] = max(best_guess[i], min_val)
    
    try:
        # Initial evaluation with error handling
        test_vals = model_func(H0, *best_guess)
        if np.any(np.isnan(test_vals)) or np.any(np.isinf(test_vals)):
            raise ValueError("Initial guess produces invalid model values")
        best_rmse = np.sqrt(np.mean((test_vals - d_delta_exp) ** 2))
    except Exception:
        best_rmse = np.inf
    
    # Define adaptive step sizes based on parameter characteristics
    step_factors = []
    for i in range(len(best_guess)):
        param_val = abs(best_guess[i])
        
        if param_val == 0:
            step_factors.append(1.0)
        elif param_val < 1:
            step_factors.append(0.1)  # Small steps for parameters < 1
        elif param_val < 100:
            step_factors.append(param_val * 0.2)  # 20% steps for medium parameters
        else:
            step_factors.append(param_val * 0.5)  # 50% steps for large parameters
    
    # Multi-scale search: coarse to fine
    scale_factors = [2.0, 1.0, 0.5, 0.2]  # Multiple scales
    
    for scale in scale_factors:
        improved_this_scale = False
        
        for iteration in range(max_iter // len(scale_factors)):
            improved_this_iter = False
            
            # Try all parameter combinations
            for i in range(len(best_guess)):
                current_step = step_factors[i] * scale
                
                # Try both directions with current step size
                for direction in [-1, 1]:
                    candidate = best_guess.copy()
                    candidate[i] += direction * current_step
                    
                    # Respect bounds with special handling for Ka parameters
                    if i < len(lower_bounds):
                        if i == 0 or i == 1:  # Ka and Kd parameters
                            min_val = max(lower_bounds[i], 1e-6)
                            candidate[i] = np.clip(candidate[i], min_val, upper_bounds[i])
                        else:
                            candidate[i] = np.clip(candidate[i], lower_bounds[i], upper_bounds[i])
                    
                    try:
                        fit_vals = model_func(H0, *candidate)
                        # Check for invalid results
                        if np.any(np.isnan(fit_vals)) or np.any(np.isinf(fit_vals)):
                            continue
                            
                        rmse = np.sqrt(np.mean((fit_vals - d_delta_exp) ** 2))
                        
                        if rmse < best_rmse:
                            best_guess = candidate.copy()
                            best_rmse = rmse
                            improved_this_iter = True
                            improved_this_scale = True
                            
                    except Exception:
                        continue
            
            # If no improvement in this iteration, try smaller steps
            if not improved_this_iter:
                for i in range(len(step_factors)):
                    step_factors[i] *= 0.5
        
        # If no improvement at this scale, move to next scale
        if not improved_this_scale:
            continue
    
    return best_guess


def _coordinate_descent_refinement(model_func, initial_guess: np.ndarray,
                                 bounds: Tuple[List[float], List[float]],
                                 H0: np.ndarray, d_delta_exp: np.ndarray) -> np.ndarray:
    """
    Fine-tune parameters using coordinate descent with golden section search.
    
    This provides final refinement by optimizing each parameter individually
    using a more sophisticated 1D optimization approach.
    """
    from scipy.optimize import minimize_scalar
    
    best_guess = initial_guess.copy()
    lower_bounds, upper_bounds = bounds
    
    # Ensure Ka parameters have minimum values
    for i in range(len(best_guess)):
        if i < len(lower_bounds):
            if i == 0 or i == 1:  # Ka and Kd parameters
                min_val = max(lower_bounds[i], 1e-6)
                best_guess[i] = max(best_guess[i], min_val)
    
    try:
        test_vals = model_func(H0, *best_guess)
        if np.any(np.isnan(test_vals)) or np.any(np.isinf(test_vals)):
            return best_guess  # Return without refinement if initial guess is bad
        best_rmse = np.sqrt(np.mean((test_vals - d_delta_exp) ** 2))
    except Exception:
        return best_guess
    
    # Refine each parameter individually
    max_refinement_iter = 3
    
    for refinement_iter in range(max_refinement_iter):
        improved = False
        
        for i in range(len(best_guess)):
            # Define objective function for parameter i
            def param_objective(x):
                candidate = best_guess.copy()
                candidate[i] = x
                try:
                    fit_vals = model_func(H0, *candidate)
                    if np.any(np.isnan(fit_vals)) or np.any(np.isinf(fit_vals)):
                        return np.inf
                    return np.sqrt(np.mean((fit_vals - d_delta_exp) ** 2))
                except Exception:
                    return np.inf
            
            # Set search bounds for this parameter
            if i < len(lower_bounds) and i < len(upper_bounds):
                if i == 0 or i == 1:  # Ka and Kd parameters
                    min_val = max(lower_bounds[i], 1e-6)
                    param_lower = max(min_val, best_guess[i] - abs(best_guess[i]) * 0.5)
                    param_upper = min(upper_bounds[i], best_guess[i] + abs(best_guess[i]) * 0.5)
                else:
                    param_lower = max(lower_bounds[i], best_guess[i] - abs(best_guess[i]))
                    param_upper = min(upper_bounds[i], best_guess[i] + abs(best_guess[i]))
                
                # Ensure valid bounds
                if param_lower >= param_upper:
                    continue
                    
                try:
                    # Use scalar optimization for this parameter
                    result = minimize_scalar(
                        param_objective,
                        bounds=(param_lower, param_upper),
                        method='bounded',
                        options={'maxiter': 20}
                    )
                    
                    if result.success and result.fun < best_rmse:
                        best_guess[i] = result.x
                        best_rmse = result.fun
                        improved = True
                        
                except Exception:
                    continue
        
        # Stop if no improvement
        if not improved:
            break
    
    return best_guess


# Enhanced fitting functions removed - the main benefit comes from smart_initial_guess
# The complex weighting/scaling approaches showed no improvement over standard methods
# Keep it simple and focus on what works: better initial parameter estimates


def enhanced_curve_fit(model_name: str, model_func, H0: np.ndarray, d_delta_exp: np.ndarray,
                      initial_guess: list, bounds: tuple, maxfev: int = 100000) -> tuple:
    """
    Enhanced curve fitting with smart initial guessing and robust error estimation.
    
    The main improvement comes from better initial parameter estimates rather than 
    complex weighting schemes. This maintains mathematical correctness while
    improving convergence and reliability.
    
    Args:
        model_name: Name of the binding model
        model_func: Model function to fit
        H0: Host concentration array
        d_delta_exp: Experimental data
        initial_guess: Initial parameter guess (should be from smart_initial_guess)
        bounds: Parameter bounds
        maxfev: Maximum function evaluations
        
    Returns:
        Tuple of (fitted_params, covariance_matrix)
    """
    
    try:
        # Use standard curve_fit with the improved initial guess
        # The main benefit comes from smart_initial_guess, not from weighting/scaling
        params, covariance = curve_fit(
            model_func, H0, d_delta_exp,
            p0=initial_guess,
            bounds=bounds,
            maxfev=maxfev
        )
        
        return params, covariance
        
    except Exception as e:
        # If fitting fails, this indicates a more serious problem
        # that won't be solved by scaling/weighting
        raise e


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
        "shapiro": (
            "Shapiro-Wilk test indicates residuals are not normally distributed. This suggests outliers, "
            "skewness, or heavy tails that may affect parameter confidence intervals."
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
        comp_stats = row.get("composite_stats", {})
        # Extract values from composite_stats if available
        spectral = comp_stats.get("spectral_ratio") if comp_stats else None
        run_ratio = comp_stats.get("run_ratio") if comp_stats else None

        # Build section headers and metric keys
        section_headers = {
            "=== Criteria ===": ["R²", "RMSE", "Weighted RMSE", "AIC", "BIC"],
            "=== Core Tests ===": [
                "Shapiro-Wilk p",
                "Ljung-Box p",
                "RESET p", 
                "Pearson corr",
                "Spearman corr",
                "Rolling R²",
            ],
            "=== Optional Tests ===": [
                "Run ratio",
                "Spectral ratio",
                "Zero-crossings", 
                "Cook's Distance",
                "Breusch-Godfrey",
            ],
        }

        table_data = [
            ["Metric", "Value", "Acceptable Range", "Possible Range", "Interpretation"]
        ]

        for section, metrics in section_headers.items():
            full_span = f"=== {section} ===".center(95)
            table_data.append([full_span, "", "", "", ""])

            for metric_name in metrics:
                if metric_name == "R²":
                    table_data.append(
                        [
                            metric_name,
                            f"{r2:.4f}",
                            "> 0.98",
                            "0 to 1",
                            interpret_diagnostic("r2", r2, 0.98, r2 > 0.98),
                        ]
                    )
                elif metric_name == "RMSE":
                    table_data.append(
                        [
                            metric_name,
                            f"{rmse:.4f}",
                            "Low",
                            "0 to ∞",
                            "Lower indicates better fit",
                        ]
                    )
                elif metric_name == "Weighted RMSE" and wrmse is not None:
                    table_data.append(
                        [
                            metric_name,
                            f"{wrmse:.4f}",
                            "Low",
                            "0 to ∞",
                            "Used for comparison across datasets",
                        ]
                    )
                elif metric_name == "AIC":
                    table_data.append(
                        [
                            metric_name,
                            f"{aic:.2f}",
                            "Low",
                            "−∞ to ∞",
                            "Used to compare model fit (penalized)",
                        ]
                    )
                elif metric_name == "BIC":
                    table_data.append(
                        [
                            metric_name,
                            f"{bic:.2f}",
                            "Low",
                            "−∞ to ∞",
                            "Stricter than AIC for model comparison",
                        ]
                    )
                elif metric_name == "Normality test":
                    shapiro_p = row.get("shapiro_p")
                    normality_pass = row.get("normality_pass", True)
                    if shapiro_p is not None:
                        table_data.append(
                            [
                                "Shapiro-Wilk p",
                                f"{shapiro_p:.3f}",
                                "> 0.05",
                                "0 to 1",
                                interpret_diagnostic(
                                    "shapiro", shapiro_p, 0.05, shapiro_p > 0.05
                                ),
                            ]
                        )
                    else:
                        table_data.append(
                            [
                                metric_name,
                                "Passed" if normality_pass else "Failed",
                                "Passed",
                                "Passed / Failed",
                                interpret_diagnostic(
                                    "normality", None, None, normality_pass
                                ),
                            ]
                        )
                elif metric_name == "Ljung-Box p" and ljung is not None:
                    table_data.append(
                        [
                            metric_name,
                            f"{ljung:.3f}",
                            "> 0.05",
                            "0 to 1",
                            interpret_diagnostic("ljung", ljung, 0.05, ljung > 0.05),
                        ]
                    )
                elif metric_name == "RESET p" and reset_p is not None:
                    table_data.append(
                        [
                            metric_name,
                            f"{reset_p:.4f}",
                            "> 0.05",
                            "0 to 1",
                            interpret_diagnostic(
                                "reset", reset_p, 0.05, reset_p > 0.05
                            ),
                        ]
                    )
                elif metric_name == "Pearson corr" and "pearson_corr" in comp_stats:
                    p = comp_stats["pearson_corr"]
                    table_data.append(
                        [
                            metric_name,
                            f"{p:.2f}",
                            "< 0.35",
                            "-1 to 1",
                            interpret_diagnostic("pearson", p, 0.35, abs(p) < 0.35),
                        ]
                    )
                elif metric_name == "Spearman corr" and "spearman_corr" in comp_stats:
                    s = comp_stats["spearman_corr"]
                    table_data.append(
                        [
                            metric_name,
                            f"{s:.2f}",
                            "< 0.35",
                            "-1 to 1",
                            interpret_diagnostic("spearman", s, 0.35, abs(s) < 0.35),
                        ]
                    )
                elif metric_name == "Rolling R²" and "avg_rolling_r2" in comp_stats:
                    r = comp_stats["avg_rolling_r2"]
                    table_data.append(
                        [
                            metric_name,
                            f"{r:.2f}",
                            "< 0.35",
                            "0 to 1",
                            interpret_diagnostic("rolling_r2", r, 0.35, r < 0.35),
                        ]
                    )
                elif metric_name == "Spectral ratio" and spectral is not None:
                    table_data.append(
                        [
                            metric_name,
                            f"{spectral:.2f}",
                            "< 0.3",
                            "0 to 1",
                            interpret_diagnostic(
                                "spectral_ratio", spectral, 0.3, spectral < 0.3
                            ),
                        ]
                    )
                elif metric_name == "Run ratio" and run_ratio is not None:
                    passed = 0.65 < run_ratio < 1.35
                    table_data.append(
                        [
                            metric_name,
                            f"{run_ratio:.2f}",
                            "≈ 1.0",
                            "0 to ~2.0",
                            interpret_diagnostic("run_ratio", run_ratio, None, passed),
                        ]
                    )
                elif metric_name == "Zero-crossings" and crossing_sim is not None:
                    passed = crossing_sim > 80
                    table_data.append(
                        [
                            metric_name,
                            f"{crossing_sim:.1f}%",
                            "> 80%",
                            "0 to 100%",
                            interpret_diagnostic(
                                "zero_crossings", crossing_sim, 80, passed
                            ),
                        ]
                    )

                elif metric_name == "Cook's Distance":
                    # This should always execute when Cook's Distance is in the Optional Tests section
                    cooks = row.get("cooks_max")
                    cooks_extreme = row.get("cooks_extreme") 
                    logging.info(f"Processing Cook's Distance: cooks_max={cooks}, cooks_extreme={cooks_extreme}")
                    
                    if cooks is not None:
                        passed = cooks < 1.0
                        display_value = f"{cooks:.3f}"
                        if cooks_extreme is not None:
                            display_value += f" ({cooks_extreme} extreme)"
                        
                        table_data.append(
                            [
                                metric_name,
                                display_value,
                                "< 1.0",
                                "0 to ∞",
                                interpret_diagnostic("cooks_distance", cooks, 1.0, passed),
                            ]
                        )
                    else:
                        table_data.append(
                            [
                                metric_name,
                                "N/A",
                                "< 1.0",
                                "0 to ∞",
                                "Cook's Distance not calculated (insufficient data points or error)",
                            ]
                        )

                elif metric_name == "Breusch-Godfrey" and bg_white is not None:
                    bg_test_name = row.get("bg_test", "Breusch-Godfrey")
                    passed = bg_white > 0.05
                    table_data.append(
                        [
                            f"{bg_test_name} p",
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

    # Shapiro-Wilk normality test (as mentioned in README)
    from scipy import stats
    try:
        if len(residuals) >= 3 and len(residuals) <= 5000:  # Shapiro-Wilk range limits
            shapiro_stat, shapiro_p = stats.shapiro(residuals)
            normality_pass = shapiro_p > 0.05
        else:
            # Fallback to 3-sigma rule for very small or very large datasets
            residuals_std = (residuals - np.mean(residuals)) / np.std(residuals)
            normality_pass = np.all(np.abs(residuals_std) < 3)
            shapiro_stat, shapiro_p = None, None
    except Exception:
        # Fallback to 3-sigma rule
        residuals_std = (residuals - np.mean(residuals)) / np.std(residuals)
        normality_pass = np.all(np.abs(residuals_std) < 3)
        shapiro_stat, shapiro_p = None, None

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
        "shapiro_stat": shapiro_stat,
        "shapiro_p": shapiro_p,
        **autocorr,
        **custom_results,
    }

    return results


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

                # Use enhanced curve fitting with weighting and scaling
                try:
                    params, cov = enhanced_curve_fit(
                        model_name, func, H0, d_delta_exp, guess, bounds, maxfev
                    )
                    fit_vals = func(H0, *params)
                    residuals = fit_vals - d_delta_exp
                    std_err = np.sqrt(np.diag(cov))
                    n_iter = maxfev  # Enhanced method doesn't track iterations directly
                    logging.info(f"Enhanced fitting completed successfully.")
                    
                except Exception as e:
                    # Fallback to least_squares if enhanced method fails
                    logging.warning(f"Enhanced fitting failed, using least_squares: {e}")
                    
                    def residuals_func(params):
                        return func(H0, *params) - d_delta_exp

                    result = least_squares(
                        residuals_func, x0=guess, bounds=bounds, max_nfev=maxfev
                    )
                    params = result.x
                    fit_vals = func(H0, *params)
                    residuals = fit_vals - d_delta_exp
                    n_iter = result.nfev * len(H0)

                    logging.info(f"Fallback fitting completed in {n_iter} function evaluations.")

                    # Estimating covariance (approximation)
                    from scipy.linalg import svd

                    _, s, VT = svd(result.jac, full_matrices=False)
                    threshold = np.finfo(float).eps * max(result.jac.shape) * s[0]
                    s = s[s > threshold]
                    VT = VT[: s.size]
                    cov = np.dot(VT.T / s**2, VT)
                    std_err = np.sqrt(np.diag(cov))

                delta_max = (
                    np.max(np.abs(d_delta_exp))
                    if np.max(np.abs(d_delta_exp)) > 0
                    else 1
                )
                normalized_residuals = residuals / global_delta_range
                weighted_rmse = np.sqrt(np.mean(normalized_residuals**2))

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

                # Create parameter dictionary with names
                param_names = model.get(
                    "parameter_names", [f"param_{i+1}" for i in range(len(params))]
                )
                named_parameters = {
                    name: value for name, value in zip(param_names, params)
                }
                named_errors = {
                    f"{name}_error": error for name, error in zip(param_names, std_err)
                }

                output_rows.append(
                    {
                        "file": filename,
                        "model": model_name,
                        "parameters": params.tolist(),
                        "parameter_names": param_names,
                        "named_parameters": named_parameters,
                        "standard_errors": std_err.tolist(),
                        "named_errors": named_errors,
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
                        **diagnostics,  # This should include cooks_max and cooks_extreme
                    }
                )
                
                # Add debug logging to verify Cook's Distance is in diagnostics
                logging.debug(f"Diagnostics keys: {list(diagnostics.keys())}")
                logging.debug(f"Cook's max: {diagnostics.get('cooks_max')}")
                logging.debug(f"Cook's extreme: {diagnostics.get('cooks_extreme')}")

                logging.info(f"Model {model_name} fit completed")
                # Log parameters with proper names
                param_names = model.get(
                    "parameter_names", [f"param_{i+1}" for i in range(len(params))]
                )
                for name, value, error in zip(param_names, params, std_err):
                    logging.info(f"  {name}: {value:.6f} ± {error:.6f}")
                logging.debug(f"Parm - results: {results}")

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
