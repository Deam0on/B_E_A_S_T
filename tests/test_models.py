import numpy as np
from scipy.optimize import curve_fit
from binding_analysis.models import (
    binding_isotherm_1_1,
    binding_isotherm_1_2,
    binding_isotherm_2_1,
    binding_dimer,
    multi_model
)

# Realistic input concentrations based on real dataset
H0 = np.linspace(0.0, 0.015, 20)
G0 = np.full_like(H0, 0.000758)

# Centralized parameter values
PARAMS = {
    "1:1": [500, 300],
    "1:2": [1e5, 1e5, 50, 200],
    "2:1": [300, 100, 50, 150],
    "dimer": [200, 100, 75, 250],
    "multi": [100, 100, 300, 20, 100, 200],
}

def test_binding_isotherm_1_1():
    Ka, d_inf = PARAMS["1:1"]
    d_delta = binding_isotherm_1_1(H0, G0, [0], Ka, d_inf)

    def model(H, Ka, d_inf):
        return binding_isotherm_1_1(H, G0, d_delta, Ka, d_inf)

    popt, _ = curve_fit(model, H0, d_delta, p0=[400, 250], maxfev=5000)
    Ka_est, d_inf_est = popt

    assert np.all(np.isfinite(popt))
    assert np.isclose(Ka_est, Ka, rtol=0.3)
    assert np.isclose(d_inf_est, d_inf, rtol=0.3)


def test_binding_isotherm_1_2():
    Ka, Kd, d1, d2 = PARAMS["1:2"]
    d_delta = binding_isotherm_1_2(H0, G0, Ka, Kd, d1, d2)

    def model(H, Ka, Kd, d1, d2):
        return binding_isotherm_1_2(H, G0, Ka, Kd, d1, d2)

    popt, _ = curve_fit(model, H0, d_delta, p0=[1e4, 1e4, 40, 180], maxfev=10000)

    assert np.all(np.isfinite(popt))


def test_binding_isotherm_2_1():
    Ka, Kd, d1, d2 = PARAMS["2:1"]
    d_delta = binding_isotherm_2_1(H0, G0, Ka, Kd, d1, d2)

    def model(H, Ka, Kd, d1, d2):
        return binding_isotherm_2_1(H, G0, Ka, Kd, d1, d2)

    popt, _ = curve_fit(model, H0, d_delta, p0=[250, 80, 40, 140], maxfev=5000)

    assert np.all(np.isfinite(popt))


def test_binding_dimer():
    Ka, Kd, d1, d2 = PARAMS["dimer"]
    d_delta = binding_dimer(H0, G0, Ka, Kd, d1, d2)

    def model(H, Ka, Kd, d1, d2):
        return binding_dimer(H, G0, Ka, Kd, d1, d2)

    popt, _ = curve_fit(model, H0, d_delta, p0=[150, 80, 50, 200], maxfev=10000)

    assert np.all(np.isfinite(popt))


def test_multi_model():
    KHG, Kd, KH2G, dG, dHG, dH2G = PARAMS["multi"]
    d_delta = multi_model(H0, G0, KHG, Kd, KH2G, dG, dHG, dH2G)

    def model(H, KHG, Kd, KH2G, dG, dHG, dH2G):
        return multi_model(H, G0, KHG, Kd, KH2G, dG, dHG, dH2G)

    popt, _ = curve_fit(model, H0, d_delta, p0=[80, 80, 200, 10, 90, 160], maxfev=10000)

    assert np.all(np.isfinite(popt))
