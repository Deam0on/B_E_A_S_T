
import os
import pandas as pd
import logging
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.diagnostic import acorr_ljungbox, acorr_breusch_godfrey, het_white

def autocorrelation_tests(H0, residuals, model_name, lags=10):
    """
    Performs Ljung-Box and Breusch-Godfrey (or White's) tests on residuals.
    Logs detailed results and interpretation.
    """
    # Ljung-Box Test
    lb_test = acorr_ljungbox(residuals, lags=[min(lags, len(H0)-2)], return_df=True)
    lb_stat = lb_test["lb_stat"].values[0]
    lb_p = lb_test["lb_pvalue"].values[0]

    # Breusch-Godfrey Test or fallback to White's
    X = sm.add_constant(H0)
    try:
        model = sm.OLS(residuals, X).fit()
        bg_stat, bg_p, _, _ = acorr_breusch_godfrey(model, nlags=min(lags, len(H0)-2))
        bg_name = "Breusch-Godfrey"
    except Exception:
        model = sm.OLS(residuals, X).fit()
        bg_stat, bg_p, _, _ = het_white(model.resid, model.model.exog)
        bg_name = "White's Test"

    ljung_fail = lb_p < 0.05
    bg_fail = bg_p < 0.05

    logging.info("-" * 70)
    logging.info(f"Autocorrelation tests for model: {model_name}")
    logging.info(f"Ljung-Box:     stat = {lb_stat:.3f}, p = {lb_p:.3f}")
    logging.info(f"{bg_name}: stat = {bg_stat:.3f}, p = {bg_p:.3f}")

    if ljung_fail and bg_fail:
        logging.warning("⚠️ Both tests detected autocorrelation! Model may need revision.")
    elif ljung_fail or bg_fail:
        logging.warning("⚠️ One test detected autocorrelation. Consider reviewing residuals.")
    else:
        logging.info("✅ No significant autocorrelation detected in residuals.")

    logging.info("-" * 70)


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
        for i, (fitted, resid, ratio) in enumerate(zip(result["fitted_values"], result["residuals"], result["H_over_G"])):
            rows.append({
                "file": result["file"],
                "model": result["model"],
                "H/G": ratio,
                "fitted_value": fitted,
                "residual": resid,
            })
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
                "RMSE": result["RMSE"]
            })

    pd.DataFrame(rows).to_csv(output_file, index=False)
