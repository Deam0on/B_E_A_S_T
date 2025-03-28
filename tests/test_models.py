import numpy as np
from scipy.optimize import curve_fit

from binding_analysis.models import (
    binding_isotherm_1_1,
    binding_isotherm_1_2,
    binding_isotherm_2_1,
    binding_dimer,
    multi_model
)

def test_binding_isotherm_1_1():
    Ka_true = 200
    d_inf_true = 250
    G0 = np.full(20, 0.01)
    H0 = np.linspace(0.001, 0.05, 20)

    d_delta = binding_isotherm_1_1(H0, G0, [0], Ka_true, d_inf_true)

    def model(H0, Ka, d_inf):
        return binding_isotherm_1_1(H0, G0, d_delta, Ka, d_inf)

    popt, _ = curve_fit(model, H0, d_delta, p0=[100, 100])
    Ka_est, d_inf_est = popt

    assert np.isclose(Ka_est, Ka_true, rtol=0.1)
    assert np.isclose(d_inf_est, d_inf_true, rtol=0.1)


def test_binding_isotherm_1_2():
    Ka_true = 100
    Kd_true = 100
    d_inf_1 = 150
    d_inf_2 = 300
    G0 = np.full(20, 0.01)
    H0 = np.linspace(0.001, 0.05, 20)

    d_delta = binding_isotherm_1_2(H0, G0, Ka_true, Kd_true, d_inf_1, d_inf_2)

    def model(H0, Ka, Kd, d1, d2):
        return binding_isotherm_1_2(H0, G0, Ka, Kd, d1, d2)

    popt, _ = curve_fit(model, H0, d_delta, p0=[50, 50, 100, 200])
    Ka_est, Kd_est, d1_est, d2_est = popt

    assert np.isclose(Ka_est, Ka_true, rtol=0.2)
    assert np.isclose(Kd_est, Kd_true, rtol=0.2)
    assert np.isclose(d1_est, d_inf_1, rtol=0.2)
    assert np.isclose(d2_est, d_inf_2, rtol=0.2)


def test_binding_isotherm_2_1():
    Ka_true = 100
    Kd_true = 100
    d_inf_1 = 100
    d_inf_2 = 200
    G0 = np.linspace(0.001, 0.05, 20)
    H0 = np.full(20, 0.01)

    d_delta = binding_isotherm_2_1(H0, G0, Ka_true, Kd_true, d_inf_1, d_inf_2)

    def model(H0, Ka, Kd, d1, d2):
        return binding_isotherm_2_1(H0, G0, Ka, Kd, d1, d2)

    popt, _ = curve_fit(model, H0, d_delta, p0=[50, 50, 50, 50])
    Ka_est, Kd_est, d1_est, d2_est = popt

    assert np.isclose(Ka_est, Ka_true, rtol=0.2)
    assert np.isclose(Kd_est, Kd_true, rtol=0.2)
    assert np.isclose(d1_est, d_inf_1, rtol=0.2)
    assert np.isclose(d2_est, d_inf_2, rtol=0.2)


def test_binding_dimer():
    Ka_true = 50
    Kd_true = 50
    d_inf_1 = 100
    d_inf_2 = 150
    G0 = np.linspace(0.001, 0.05, 20)
    H0 = np.full(20, 0.01)

    d_delta = binding_dimer(H0, G0, Ka_true, Kd_true, d_inf_1, d_inf_2)

    def model(H0, Ka, Kd, d1, d2):
        return binding_dimer(H0, G0, Ka, Kd, d1, d2)

    popt, _ = curve_fit(model, H0, d_delta, p0=[25, 25, 75, 75])
    Ka_est, Kd_est, d1_est, d2_est = popt

    assert np.isclose(Ka_est, Ka_true, rtol=0.2)
    assert np.isclose(Kd_est, Kd_true, rtol=0.2)
    assert np.isclose(d1_est, d_inf_1, rtol=0.2)
    assert np.isclose(d2_est, d_inf_2, rtol=0.2)


def test_multi_model():
    KHG = 100
    Kd = 100
    KH2G = 200
    dG = 20
    dHG = 100
    dH2G = 180
    G0 = np.full(20, 0.01)
    H0 = np.linspace(0.001, 0.05, 20)

    d_delta = multi_model(H0, G0, KHG, Kd, KH2G, dG, dHG, dH2G)

    def model(H0, KHG, Kd, KH2G, dG, dHG, dH2G):
        return multi_model(H0, G0, KHG, Kd, KH2G, dG, dHG, dH2G)

    popt, _ = curve_fit(model, H0, d_delta, p0=[80, 80, 180, 10, 90, 160])
    KHG_e, Kd_e, KH2G_e, dG_e, dHG_e, dH2G_e = popt

    assert np.isclose(KHG_e, KHG, rtol=0.2)
    assert np.isclose(Kd_e, Kd, rtol=0.2)
    assert np.isclose(KH2G_e, KH2G, rtol=0.2)
    assert np.isclose(dG_e, dG, rtol=0.2)
    assert np.isclose(dHG_e, dHG, rtol=0.2)
    assert np.isclose(dH2G_e, dH2G, rtol=0.2)
