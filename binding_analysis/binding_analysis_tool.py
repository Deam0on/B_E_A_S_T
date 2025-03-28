import sys
import os
import logging
from analysis import process_csv_files_in_folder
from utils import delete_old_result_files
from datetime import datetime

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

def main():
    setup_logging()
    input_folder = "data_input"
    output_folder = "results"

    os.makedirs(output_folder, exist_ok=True)
    delete_old_result_files(output_folder)
    process_csv_files_in_folder(input_folder, output_folder)

    logging.info("Analysis completed. Results and log saved to /results.")

if __name__ == "__main__":
    main()
