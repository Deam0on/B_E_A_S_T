
import os
import pandas as pd
import logging

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
