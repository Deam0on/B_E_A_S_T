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
    "1:1": [1e4, 100, 350],
    "1:2": [1e4, 1e3, 150, 350],
    "2:1": [1e4, 1e3, 150, 350],
    "dimer": [1e3, 1e3, 150, 350],
    "multi": [1e3, 1e3, 1e3, 150, 300, 400],
}

def validate_fit(model_lambda, H0, true_params, initial_guess, bounds, rtol=1e-2):
    d_calc = model_lambda(H0, *true_params)
    popt, _ = curve_fit(model_lambda, H0, d_calc, p0=initial_guess, bounds=bounds, maxfev=100000)
    fit = model_lambda(H0, *popt)

    if not np.allclose(fit, d_calc, rtol=rtol):
        print("Expected:", d_calc)
        print("Obtained:", fit)
        print("Params:", popt)
        return False
    return True

def test_model_1_1():
    model = model_definitions(H0, G0, d_delta_exp)["1:1"]
    assert validate_fit(model["lambda"], H0, PARAMS["1:1"], model["initial_guess"], model["bounds"])

def test_model_1_2():
    model = model_definitions(H0, G0, d_delta_exp)["1:2"]
    assert validate_fit(model["lambda"], H0, PARAMS["1:2"], model["initial_guess"], model["bounds"])

def test_model_2_1():
    model = model_definitions(H0, G0, d_delta_exp)["2:1"]
    assert validate_fit(model["lambda"], H0, PARAMS["2:1"], model["initial_guess"], model["bounds"])

def test_model_dimer():
    model = model_definitions(H0, G0, d_delta_exp)["dimer"]
    assert validate_fit(model["lambda"], H0, PARAMS["dimer"], model["initial_guess"], model["bounds"])

def test_model_multi():
    model = model_definitions(H0, G0, d_delta_exp)["multi"]
    assert validate_fit(model["lambda"], H0, PARAMS["multi"], model["initial_guess"], model["bounds"])
