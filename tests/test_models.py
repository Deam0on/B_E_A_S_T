import numpy as np
from scipy.optimize import curve_fit
from binding_analysis.models import model_definitions

# Shared test input
H0 = np.linspace(0, 0.015, 20)
G0 = np.full_like(H0, 0.01)
d_delta_exp = np.linspace(0, 2.5e3, len(H0))

# Expected realistic parameter guesses for good convergence
PARAMS = {
    "1:1": [1e4, 300],
    "1:2": [1e4, 1e3, 150, 300],
    "2:1": [1e4, 1e3, 150, 300],
    "dimer": [1e4, 1e3, 150, 300],
    "multi": [1e4, 1e3, 1e4, 2700, 150, 350]
}

def test_model_1_1():
    model = model_definitions(H0, G0)["1:1"]
    d_sim = model["lambda"](H0, *PARAMS["1:1"])
    popt, _ = curve_fit(model["lambda"], H0, d_sim, p0=model["initial_guess"], bounds=model["bounds"])
    assert len(popt) == 2

def test_model_1_2():
    model = model_definitions(H0, G0)["1:2"]
    d_sim = model["lambda"](H0, *PARAMS["1:2"])
    popt, _ = curve_fit(model["lambda"], H0, d_sim, p0=model["initial_guess"], bounds=model["bounds"])
    assert len(popt) == 4

def test_model_2_1():
    model = model_definitions(H0, G0)["2:1"]
    d_sim = model["lambda"](H0, *PARAMS["2:1"])
    popt, _ = curve_fit(model["lambda"], H0, d_sim, p0=model["initial_guess"], bounds=model["bounds"])
    assert len(popt) == 4

def test_model_dimer():
    model = model_definitions(H0, G0)["dimer"]
    d_sim = model["lambda"](H0, *PARAMS["dimer"])
    popt, _ = curve_fit(model["lambda"], H0, d_sim, p0=model["initial_guess"], bounds=model["bounds"])
    assert len(popt) == 4

def test_model_multi():
    model = model_definitions(H0, G0)["multi"]
    d_sim = model["lambda"](H0, *PARAMS["multi"])
    popt, _ = curve_fit(model["lambda"], H0, d_sim, p0=model["initial_guess"], bounds=model["bounds"])
    assert len(popt) == 6
