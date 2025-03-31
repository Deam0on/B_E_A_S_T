import numpy as np
from scipy.optimize import curve_fit
from binding_analysis.models import (
    binding_isotherm_1_1,
    binding_isotherm_1_2,
    binding_isotherm_2_1,
    binding_dimer,
    multi_model
)

# Realistic input concentrations
H0 = np.linspace(0.0, 0.015, 20)
G0 = np.full_like(H0, 0.000758)

# Centralized test parameters
PARAMS = {
    "1:1": [500, 300, 2800],                     # Ka, dG, dHG
    "1:2": [1e5, 1e5, 2800, 100, 400],           # K1, K2, dG, dHG, dHG2
    "2:1": [1e5, 1e5, 2800, 100, 400],           # K1, K2, dG, dHG, dH2G
    "dimer": [200, 100, 2800, 350],              # Ka, Kd, dG, dHG
    "multi": [100, 100, 300, 2800, 100, 400]     # KHG, Kd, KH2G, dG, dHG, dH2G
}

def test_binding_isotherm_1_1():
    Ka, dG, dHG = PARAMS["1:1"]
    d_delta = binding_isotherm_1_1(H0, G0, Ka, dG, dHG)

    def model(H, Ka, dG, dHG):
        return binding_isotherm_1_1(H, G0, Ka, dG, dHG)

    popt, _ = curve_fit(model, H0, d_delta, p0=[400, 2700, 3000], maxfev=5000)
    assert np.all(np.isfinite(popt))


def test_binding_isotherm_1_2():
    K1, K2, dG, dHG, dHG2 = PARAMS["1:2"]
    d_delta = binding_isotherm_1_2(H0, G0, K1, K2, dG, dHG, dHG2)

    def model(H, K1, K2, dG, dHG, dHG2):
        return binding_isotherm_1_2(H, G0, K1, K2, dG, dHG, dHG2)

    popt, _ = curve_fit(model, H0, d_delta, p0=[1e4, 1e4, 2700, 150, 350], maxfev=10000)
    assert np.all(np.isfinite(popt))


def test_binding_isotherm_2_1():
    K1, K2, dG, dHG, dH2G = PARAMS["2:1"]
    d_delta = binding_isotherm_2_1(H0, G0, K1, K2, dG, dHG, dH2G)

    def model(H, K1, K2, dG, dHG, dH2G):
        return binding_isotherm_2_1(H, G0, K1, K2, dG, dHG, dH2G)

    popt, _ = curve_fit(model, H0, d_delta, p0=[1e4, 1e4, 2700, 150, 350], maxfev=10000)
    assert np.all(np.isfinite(popt))


def test_binding_dimer():
    Ka, Kd, dG, dHG = PARAMS["dimer"]
    d_delta = binding_dimer(H0, G0, Ka, Kd, dG, dHG)

    def model(H, Ka, Kd, dG, dHG):
        return binding_dimer(H, G0, Ka, Kd, dG, dHG)

    popt, _ = curve_fit(model, H0, d_delta, p0=[150, 80, 2700, 300], maxfev=10000)
    assert np.all(np.isfinite(popt))


def test_multi_model():
    KHG, Kd, KH2G, dG, dHG, dH2G = PARAMS["multi"]
    d_delta = multi_model(H0, G0, KHG, Kd, KH2G, dG, dHG, dH2G)

    def model(H, KHG, Kd, KH2G, dG, dHG, dH2G):
        return multi_model(H, G0, KHG, Kd, KH2G, dG, dHG, dH2G)

    popt, _ = curve_fit(model, H0, d_delta, p0=[80, 80, 200, 2700, 100, 350], maxfev=15000)
    assert np.all(np.isfinite(popt))
