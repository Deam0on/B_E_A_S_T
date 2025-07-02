"""
BEAST - Binding Evaluation and Analysis Software Tool

This module provides a command-line interface for analyzing NMR titration data
and fitting it to various thermodynamic binding models.

Author: Filip Hládek
License: MIT
"""

import argparse
import logging
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict

import yaml
from analysis import process_csv_files_in_folder
from utils import delete_old_result_files


def load_config(filename: str = "config.yaml") -> Dict[str, Any]:
    """
    Load configuration from a YAML file.

    Args:
        filename: Path to the configuration file (default: "config.yaml")

    Returns:
        Configuration dictionary loaded from YAML file

    Raises:
        FileNotFoundError: If config file doesn't exist
        yaml.YAMLError: If config file is malformed
    """
    base_path = Path(__file__).parent
    config_path = base_path / filename

    try:
        with open(config_path, "r", encoding="utf-8") as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        logging.error(f"Configuration file '{config_path}' not found")
        raise
    except yaml.YAMLError as e:
        logging.error(f"Error parsing configuration file: {e}")
        raise


def log_uncaught_exceptions(exc_type, exc_value, exc_traceback) -> None:
    """
    Log uncaught exceptions to the logging system.

    Args:
        exc_type: Exception type
        exc_value: Exception value
        exc_traceback: Exception traceback
    """
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logging.critical(
        "Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback)
    )


sys.excepthook = log_uncaught_exceptions


def setup_logging() -> None:
    """
    Configure logging to write to both file and console.

    Creates a results directory if it doesn't exist and sets up
    logging handlers for both file output and console output.
    """
    results_dir = Path("results")
    results_dir.mkdir(exist_ok=True)
    log_path = results_dir / "log.txt"

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_path, mode="w", encoding="utf-8"),
            logging.StreamHandler(),
        ],
    )

    logging.info("Binding analysis started")
    logging.info(f"Timestamp: {datetime.now().isoformat()}")


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments.

    Returns:
        Parsed command line arguments
    """
    parser = argparse.ArgumentParser(
        description="BEAST - Binding Evaluation and Analysis Software Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                                    # Use default config
  %(prog)s --input_dir data --output_dir out  # Custom directories
  %(prog)s --skip_tests --no_normalized       # Disable tests and normalized plots
        """,
    )
    parser.add_argument(
        "--config",
        type=str,
        default="config.yaml",
        help="Path to config YAML file (default: %(default)s)",
    )
    parser.add_argument(
        "--skip_tests", action="store_true", help="Disable residual diagnostic tests"
    )
    parser.add_argument(
        "--no_normalized",
        action="store_true",
        help="Do not show/save normalized residual plots",
    )
    parser.add_argument(
        "--input_dir", type=str, help="Override input directory (default: from config)"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        help="Override output directory (default: from config)",
    )
    parser.add_argument(
        "--global_norm",
        action="store_true",
        help="Normalize residuals using global max Δδ",
    )
    parser.add_argument(
        "--custom_residual_check",
        action="store_true",
        help="Enable custom pattern-based residual test",
    )
    return parser.parse_args()


def main() -> None:
    """
    Main entry point for the binding analysis tool.

    Sets up logging, parses arguments, loads configuration, and runs the analysis.
    """
    setup_logging()
    args = parse_args()

    try:
        config = load_config(args.config)
    except (FileNotFoundError, yaml.YAMLError) as e:
        logging.error(f"Failed to load configuration: {e}")
        sys.exit(1)

    # Override input/output directory if provided via CLI
    input_dir = args.input_dir or config.get("general", {}).get(
        "input_dir", "data_input"
    )
    output_dir = args.output_dir or config.get("general", {}).get(
        "results_dir", "results"
    )

    # Validate input directory exists
    input_path = Path(input_dir)
    if not input_path.exists():
        logging.error(f"Input directory '{input_path}' does not exist")
        sys.exit(1)

    # Update config with final values
    config["general"]["input_dir"] = str(input_path)
    config["general"]["results_dir"] = output_dir

    # Store CLI flags in config
    config["cli_flags"] = {
        "skip_tests": args.skip_tests,
        "no_normalized": args.no_normalized,
        "custom_residual_check": args.custom_residual_check,
    }

    # Create output directory and clean old results
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    delete_old_result_files(str(output_path))

    # Run the analysis
    try:
        process_csv_files_in_folder(
            config, skip_tests=args.skip_tests, plot_normalized=not args.no_normalized
        )

        logging.info(
            f"Analysis completed successfully. Results saved to '{output_path}'."
        )

    except Exception as e:
        logging.error(f"Analysis failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
