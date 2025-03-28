
import os
import logging
from analysis import process_csv_files_in_folder
from utils import delete_old_result_files

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='[%(levelname)s] %(message)s'
    )

def main():
    setup_logging()
    input_folder = "data_input"
    output_folder = "results"

    os.makedirs(output_folder, exist_ok=True)
    delete_old_result_files(output_folder)
    process_csv_files_in_folder(input_folder, output_folder)

if __name__ == "__main__":
    main()
