"""
Setup script for BEAST - Binding Evaluation and Analysis Software Tool
"""

from pathlib import Path

from setuptools import find_packages, setup

# Read the README file
this_directory = Path(__file__).parent
long_description = (
    (this_directory / "README.md").read_text(encoding="utf-8")
    if (this_directory / "README.md").exists()
    else ""
)

# Read requirements
requirements = []
req_file = this_directory / "binding_analysis" / "requirements.txt"
if req_file.exists():
    requirements = req_file.read_text().strip().split("\n")
    requirements = [
        req.strip() for req in requirements if req.strip() and not req.startswith("#")
    ]

setup(
    name="beast-binding-analysis",
    version="1.1.9",
    author="Filip Hládek",
    author_email="info.f@hladek.cz",
    license="MIT",
    description="BEAST - Binding Evaluation and Analysis Software Tool for NMR titration data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Deam0on/B_E_A_S_T",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Information Analysis",
    ],
    python_requires=">=3.7",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.12",
            "black>=21.0",
            "flake8>=3.9",
            "mypy>=0.910",
        ],
    },
    entry_points={
        "console_scripts": [
            "beast=binding_analysis.binding_analysis_tool:main",
        ],
    },
    include_package_data=True,
    package_data={
        "binding_analysis": ["*.yaml", "*.yml"],
    },
    keywords="nmr titration binding thermodynamics curve-fitting chemistry",
    project_urls={
        "Bug Reports": "https://github.com/Deam0on/B_E_A_S_T/issues",
        "Source": "https://github.com/Deam0on/B_E_A_S_T",
        "Documentation": "https://github.com/Deam0on/B_E_A_S_T/wiki",
    },
)
