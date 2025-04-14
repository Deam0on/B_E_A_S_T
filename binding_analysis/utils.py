import os
import pandas as pd
import logging
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.diagnostic import acorr_ljungbox, acorr_breusch_godfrey, het_white
from scipy.optimize import curve_fit
from scipy.stats import skew, kurtosis, normaltest
from scipy.stats import pearsonr, spearmanr
from statsmodels.sandbox.stats.runs import runstest_1samp
from scipy.fft import rfft

def custom_residual_pattern_test(residuals):
    if len(residuals) < 9:
        return {"composite_stats": None, "composite_flagged": None}

    residuals = np.asarray(residuals)
    n = len(residuals)

    # Lag-1 Pearson and Spearman
    pearson_corr, _ = pearsonr(residuals[:-1], residuals[1:])
    spearman_corr, _ = spearmanr(residuals[:-1], residuals[1:])

    # Spectral energy ratio (low freq)
    freqs = rfft(residuals)
    energy = np.abs(freqs)**2
    low_freq_energy = np.sum(energy[:n//10])
    total_energy = np.sum(energy)
    spectral_ratio = low_freq_energy / total_energy if total_energy > 0 else 0

    # Rolling R² (linear fit on residual windows)
    window = min(5, len(residuals) - 2)
    rolling_r2 = []
    for i in range(n - window + 1):
        x = np.arange(window)
        y = residuals[i:i+window]
        slope, intercept = np.polyfit(x, y, 1)
        y_pred = slope * x + intercept
        ss_res = np.sum((y - y_pred) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        rolling_r2.append(r2)
    avg_r2 = np.mean(rolling_r2)

    # Runs test (sign changes)
    signs = np.sign(residuals - np.mean(residuals))
    run_stat, p_run = runstest_1samp(signs)
    expected_runs = (2 * np.sum(signs != 0) - 1) / 3
    actual_runs = np.sum(np.diff(signs) != 0)
    run_ratio = actual_runs / expected_runs if expected_runs else 1

    # Heuristic thresholds
    flagged = (
        abs(pearson_corr) > 0.35 or
        abs(spearman_corr) > 0.35 or
        spectral_ratio > 0.3 or
        avg_r2 > 0.35 or
        run_ratio < 0.65 or run_ratio > 1.35
    )

    return {
        "composite_stats": {
            "pearson_corr": pearson_corr,
            "spearman_corr": spearman_corr,
            "spectral_ratio": spectral_ratio,
            "avg_rolling_r2": avg_r2,
            "run_ratio": run_ratio
        },
        "composite_flagged": flagged
    }

def collect_global_max_deltadelta(input_folder: str) -> float:
    """
    Scans all CSV files in the input folder and returns the maximum Δδ across all files.
    """
    import pandas as pd
    import numpy as np
    from pathlib import Path

    max_ddeltas = []

    for csv_file in Path(input_folder).glob("*.csv"):
        try:
            df = pd.read_csv(csv_file)
            delta = np.abs(df["delta"] - df["delta"].iloc[0])
            delta.iloc[0] = 0
            max_ddeltas.append(delta.max())
        except Exception as e:
            logging.warning(f"Failed to parse {csv_file.name}: {e}")

    return max(max_ddeltas) if max_ddeltas else 1.0


def smart_initial_guess(model_func, H0, d_delta_exp, default_guess, bounds, max_iter=5, perturbation=0.1):
    best_guess = default_guess
    best_score = np.inf

    for _ in range(max_iter):
        guess = [p * (1 + np.random.uniform(-perturbation, perturbation)) for p in default_guess]
        try:
            popt, _ = curve_fit(model_func, H0, d_delta_exp, p0=guess, bounds=bounds, maxfev=100000)
            residuals = model_func(H0, *popt) - d_delta_exp
            score = np.sum(residuals ** 2)
            if score < best_score:
                best_score = score
                best_guess = guess
        except Exception:
            continue

    return best_guess


def compare_models_by_metric(output_rows, metric="AIC"):
    try:
        model_metrics = [(row["model"], row[metric]) for row in output_rows]
        sorted_models = sorted(model_metrics, key=lambda x: x[1])
        logging.info("Model Ranking by %s:", metric)
        for rank, (model, val) in enumerate(sorted_models, start=1):
            logging.info(f"{rank}. {model}: {metric} = {val:.2f}")
    except Exception as e:
        logging.warning(f"Model comparison skipped due to: {e}")


def advanced_residual_diagnostics(residuals, model_name):
    try:
        s = skew(residuals)
        k = kurtosis(residuals)
        stat, p = normaltest(residuals)
        # logging.info(f"[{model_name}] Residual Diagnostics -> Skew: {s:.2f}, Kurtosis: {k:.2f}, Normality p = {p:.4f}")
        return {"skew": s, "kurtosis": k, "normality_p": p}
    except Exception as e:
        logging.warning(f"Residual diagnostics for {model_name} skipped due to: {e}")
        return {}


def validate_data(df):
    required_cols = ['H', 'G', 'delta']
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        logging.error(f"Missing columns: {missing}")
        return False
    if df[required_cols].isnull().any().any():
        logging.warning("Detected NaNs in input data.")
    return True


def autocorrelation_tests(H0, residuals, model_name, lags=10):
    results = {}
    n = len(H0)
    X = sm.add_constant(H0)

    lags = min(10, n // 5)

    # logging.info("-" * 70)
    # logging.info(f"Residual diagnostics for model: {model_name}")

    # Ljung-Box
    if n >= 2 * lags + 1:
        lb_test = acorr_ljungbox(residuals, lags=[min(lags, n - 2)], return_df=True)
        lb_stat = lb_test["lb_stat"].values[0]
        lb_p = lb_test["lb_pvalue"].values[0]
        results.update({
            "ljung_stat": lb_stat,
            "ljung_p": lb_p,
            "ljung_failed": lb_p < 0.05
        })
        # logging.info(f"Ljung-Box:     stat = {lb_stat:.3f}, p = {lb_p:.4f}")
    else:
        logging.info("Ljung-Box test skipped (too few data points).")
        results.update({
            "ljung_stat": None,
            "ljung_p": None,
            "ljung_failed": None
        })

    # Breusch-Godfrey or White
    try:
        model = sm.OLS(residuals, X).fit()
        bg_stat, bg_p, _, _ = acorr_breusch_godfrey(model, nlags=min(lags, n - 2))
        bg_name = "Breusch-Godfrey"
    except Exception:
        model = sm.OLS(residuals, X).fit()
        bg_stat, bg_p, _, _ = het_white(model.resid, model.model.exog)
        bg_name = "White's Test"

    results.update({
        "bg_test": bg_name,
        "bg_stat": bg_stat,
        "bg_p": bg_p,
        "bg_failed": bg_p < 0.05
    })
    # logging.info(f"{bg_name}: stat = {bg_stat:.3f}, p = {bg_p:.4f}")

    # Ramsey RESET test
    if n >= 15:
        try:
            model = sm.OLS(residuals, X).fit()
            from statsmodels.stats.diagnostic import linear_reset
            reset_test = linear_reset(model, power=2, use_f=True)
            results.update({
                "reset_stat": reset_test.fvalue,
                "reset_p": reset_test.pvalue,
                "reset_failed": reset_test.pvalue < 0.05
            })
            # logging.info(f"Ramsey RESET:  stat = {reset_test.fvalue:.3f}, p = {reset_test.pvalue:.4f}")
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
            results.update({
                "cooks_max": max_cook,
                "cooks_extreme": n_extreme
            })
            # logging.info(f"Cook’s Distance: max = {max_cook:.4f}, extreme (>4/n): {n_extreme}")
        except Exception as e:
            logging.warning(f"Cook’s Distance failed: {e}")
    else:
        logging.info("Cook’s Distance skipped (too few data points).")

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
        # logging.info(f"Zero-crossings: {zero_crossings} | Similarity to white noise: {similarity:.2f}%")
    except Exception as e:
        logging.warning(f"Zero-crossing similarity test failed: {e}")

    logging.info("-" * 70)
    return results


def delete_old_result_files(folder_path):
    for f in os.listdir(folder_path):
        if f.endswith('_results.csv') or f.endswith('_plot.png'):
            try:
                os.remove(os.path.join(folder_path, f))
                logging.info(f"Deleted: {f}")
            except Exception as e:
                logging.warning(f"Failed to delete {f}: {e}")


def save_combined_csv(results, output_file):
    rows = []
    for result in results:
        num_points = len(result["fitted_values"])
        for i in range(num_points):
            row = {
                "file": result["file"],
                "model": result["model"],
                "H/G": result["H_over_G"][i],
                "fitted_value": result["fitted_values"][i],
                "residual": result["residuals"][i]
            }
            if result.get("normalized_residuals"):
                row["normalized_residual"] = result["normalized_residuals"][i]
            rows.append(row)

        for i, param in enumerate(result["parameters"]):
            rows.append({
                "file": result["file"],
                "model": result["model"],
                "parameter_index": i + 1,
                "parameter_value": param,
                "standard_error": result["standard_errors"][i],
                "confidence_interval_lower": result["confidence_intervals"][i][0],
                "confidence_interval_upper": result["confidence_intervals"][i][1],
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
                "composite_stats": str(result.get("composite_stats"))
            })

    pd.DataFrame(rows).to_csv(output_file, index=False)
