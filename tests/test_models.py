import numpy as np
from scipy.optimize import curve_fit
from binding_analysis.models import model_definitions

# Simulated test data
H0 = np.linspace(0, 0.015, 20)
G0 = np.full_like(H0, 0.01)
delta = 2800 - 2000 * (H0 / (H0 + 0.005))
d_delta_exp = np.abs(delta - delta[0])
d_delta_exp[0] = 0

PARAMS = {
    "1:1": [1e4, 350],
    "1:2": [1e4, 1e3, 150, 350],
    "2:1": [1e4, 1e3, 150, 350],
    "dimer": [1e3, 1e3, 150, 350],
    "multi": [1e3, 1e3, 1e3, 150, 300, 400],
}

def test_model_1_1():
    model = model_definitions(H0, G0, d_delta_exp)["1:1"]
    d_calc = model["lambda"](H0, *PARAMS["1:1"])
    popt, _ = curve_fit(model["lambda"], H0, d_calc, p0=model["initial_guess"], bounds=model["bounds"])
    assert np.allclose(model["lambda"](H0, *popt), d_calc, rtol=1e-2)

def test_model_1_2():
    model = model_definitions(H0, G0, d_delta_exp)["1:2"]
    d_calc = model["lambda"](H0, *PARAMS["1:2"])
    popt, _ = curve_fit(model["lambda"], H0, d_calc, p0=model["initial_guess"], bounds=model["bounds"], maxfev=100000)
    assert np.allclose(model["lambda"](H0, *popt), d_calc, rtol=1e-2)

def test_model_2_1():
    model = model_definitions(H0, G0, d_delta_exp)["2:1"]
    d_calc = model["lambda"](H0, *PARAMS["2:1"])
    popt, _ = curve_fit(model["lambda"], H0, d_calc, p0=model["initial_guess"], bounds=model["bounds"], maxfev=100000)
    assert np.allclose(model["lambda"](H0, *popt), d_calc, rtol=1e-2)

def test_model_dimer():
    model = model_definitions(H0, G0, d_delta_exp)["dimer"]
    d_calc = model["lambda"](H0, *PARAMS["dimer"])
    popt, _ = curve_fit(model["lambda"], H0, d_calc, p0=model["initial_guess"], bounds=model["bounds"], maxfev=100000)
    assert np.allclose(model["lambda"](H0, *popt), d_calc, rtol=1e-2)

def test_model_multi():
    model = model_definitions(H0, G0, d_delta_exp)["multi"]
    d_calc = model["lambda"](H0, *PARAMS["multi"])
    popt, _ = curve_fit(model["lambda"], H0, d_calc, p0=model["initial_guess"], bounds=model["bounds"], maxfev=100000)
    assert np.allclose(model["lambda"](H0, *popt), d_calc, rtol=1e-2)
