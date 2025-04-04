""import sys
import os
import logging
import argparse
from analysis import process_csv_files_in_folder
from utils import delete_old_result_files
from datetime import datetime
import yaml

def load_config(filename="config.yaml"):
    base_path = os.path.dirname(os.path.abspath(__file__))
    config_path = os.path.join(base_path, filename)
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def log_uncaught_exceptions(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logging.critical("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))

sys.excepthook = log_uncaught_exceptions

def setup_logging():
    os.makedirs("results", exist_ok=True)
    log_path = os.path.join("results", "log.txt")

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_path, mode='w', encoding='utf-8'),
            logging.StreamHandler()
        ]
    )

    logging.info("Binding analysis started")
    logging.info(f"Timestamp: {datetime.now().isoformat()}")

def parse_args():
    parser = argparse.ArgumentParser(description="Run NMR binding model analysis.")
    parser.add_argument("--config", type=str, default="config.yaml", help="Path to YAML configuration file")
    parser.add_argument("--skip_tests", action="store_true", help="Skip residual statistical diagnostics")
    parser.add_argument("--no_normalized", action="store_true", help="Skip normalized residual plotting")
    return parser.parse_args()

def main():
    setup_logging()
    args = parse_args()
    config = load_config(args.config)

    # Inject CLI flags into config
    config["cli_flags"] = {
        "skip_tests": args.skip_tests,
        "no_normalized": args.no_normalized
    }

    input_folder = config.get("general", {}).get("input_dir", "data_input")
    output_folder = config.get("general", {}).get("results_dir", "results")

    os.makedirs(output_folder, exist_ok=True)
    delete_old_result_files(output_folder)
    process_csv_files_in_folder(config)

    logging.info("Analysis completed. Results and log saved to /results.")

if __name__ == "__main__":
    main()
