import sys
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
    parser = argparse.ArgumentParser(description="Binding Isotherm CLI Tool")
    parser.add_argument("--config", type=str, default="config.yaml", help="Path to config YAML file.")
    parser.add_argument("--skip_tests", action="store_true", help="Disable residual diagnostic tests.")
    parser.add_argument("--no_normalized", action="store_true", help="Do not show/save normalized residual plots.")
    parser.add_argument("--input_dir", type=str, help="Override input directory (default: from config)")
    parser.add_argument("--output_dir", type=str, help="Override output directory (default: from config)")
    parser.add_argument("--global_norm", action="store_true", help="Normalize residuals using global max Δδ.")
    parser.add_argument("--custom_residual_check", action="store_true", help="Enable custom pattern-based residual test.")
    return parser.parse_args()

def main():
    setup_logging()
    args = parse_args()
    config = load_config(args.config)

    # Override input/output directory if provided via CLI
    input_dir = args.input_dir or config.get("general", {}).get("input_dir", "data_input")
    output_dir = args.output_dir or config.get("general", {}).get("results_dir", "results")

    # Update config with final values
    config["general"]["input_dir"] = input_dir
    config["general"]["results_dir"] = output_dir

    # Store CLI flags in config
    config["cli_flags"] = {
        "skip_tests": args.skip_tests,
        "no_normalized": args.no_normalized
    }

    os.makedirs(output_dir, exist_ok=True)
    delete_old_result_files(output_dir)

    # Pass CLI flags explicitly to the analysis function
    process_csv_files_in_folder(
        config,
        skip_tests=args.skip_tests,
        plot_normalized=not args.no_normalized
    )

    logging.info("Analysis completed. Results and log saved to /results.")


if __name__ == "__main__":
    main()
