"""
BEAST - Binding Evaluation and Analysis Software Tool

A comprehensive Python package for analyzing NMR titration data and fitting
it to various thermodynamic binding models.

Author: Filip Hládek
License: MIT
Version: 1.1.9
"""

__version__ = "1.1.9"
__author__ = "Filip Hládek"
__email__ = "info.f@hladek.cz"
__license__ = "MIT"

# Import main modules for easy access
from . import analysis, models, utils

# Import main functions
from .analysis import plot_results, process_csv_files_in_folder
from .models import model_definitions
from .utils import collect_global_max_deltadelta, validate_data

__all__ = [
    "models",
    "utils",
    "analysis",
    "process_csv_files_in_folder",
    "plot_results",
    "model_definitions",
    "validate_data",
    "collect_global_max_deltadelta",
]
