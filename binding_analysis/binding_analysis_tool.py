
import os
import logging
from analysis import process_csv_files_in_folder
from utils import delete_old_result_files

def setup_logging():
    os.makedirs("results", exist_ok=True)
    log_path = os.path.join("results", "log.txt")

    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        handlers=[
            logging.FileHandler(log_path, mode='w'),
            logging.StreamHandler()  # still print to Colab terminal
        ]
    )

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
