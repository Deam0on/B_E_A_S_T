
import numpy as np

def binding_isotherm_1_1(H0, G0, Ka, dG, dHG):
    d_obs = np.zeros_like(H0)
    epsilon = 1e-10
    for i in range(len(H0)):
        a = Ka
        b = -Ka * (H0[i] + G0[i])
        c = Ka * H0[i] * G0[i]
        roots = np.roots([a, b, c])
        real_roots = roots[np.isreal(roots)].real
        HG = np.max(real_roots[real_roots > epsilon]) if len(real_roots[real_roots > epsilon]) > 0 else epsilon
        G = G0[i] - HG
        d_obs[i] = (G * dG + HG * dHG) / G0[i]
    return d_obs

def binding_isotherm_1_2(H0, G0, K1, K2, dG, dHG, dHG2):
    d_obs = np.zeros_like(H0)
    epsilon = 1e-10
    for i in range(len(H0)):
        a = 2 * K1 * K2 * H0[i]
        b = K1 * H0[i]
        c = 1
        d = -G0[i]
        roots = np.roots([a, b, c, d])
        real_roots = roots[np.isreal(roots)].real
        G = np.min(real_roots[real_roots > epsilon]) if len(real_roots[real_roots > epsilon]) > 0 else epsilon
        HG = K1 * H0[i] * G
        HG2 = K1 * K2 * H0[i] * G**2
        free_G = G0[i] - HG - 2 * HG2
        d_obs[i] = (free_G * dG + HG * dHG + HG2 * dHG2) / G0[i]
    return d_obs

def binding_isotherm_2_1(H0, G0, K1, K2, dG, dHG, dH2G):
    d_obs = np.zeros_like(H0)
    epsilon = 1e-10
    for i in range(len(H0)):
        a = 2 * K1 * K2 * G0[i]
        b = K1 * G0[i]
        c = 1
        d = -H0[i]
        roots = np.roots([a, b, c, d])
        real_roots = roots[np.isreal(roots)].real
        H = np.min(real_roots[real_roots > epsilon]) if len(real_roots[real_roots > epsilon]) > 0 else epsilon
        HG = K1 * H * G0[i]
        H2G = K1 * K2 * H**2 * G0[i]
        G = G0[i] - HG - H2G
        d_obs[i] = (G * dG + HG * dHG + H2G * dH2G) / G0[i]
    return d_obs


def binding_dimer(H0, G0, Ka, Kd, d_inf_1, d_inf_2):
    d_obs = np.zeros_like(H0)
    epsilon = 1e-10
    for i in range(len(H0)):
        a = -2 * Kd
        b = -1 * (G0[i] * Ka + 1)
        c = H0[i]
        roots = np.roots([a, b, c])
        real_roots = roots[np.isreal(roots)].real
        H_free = np.min(real_roots[real_roots > epsilon]) if len(real_roots[real_roots > epsilon]) > 0 else epsilon
        numerator = (d_inf_1 * G0[i]) / (1 + Ka * H_free) + (d_inf_2 * Ka * H_free * G0[i])
        denominator = (G0[i] / (1 + Ka * H_free)) + (Ka * H_free * G0[i])
        d_obs[i] = numerator / denominator
    return d_obs

def multi_model(H0, G0, KHG, Kd, KH2G, dG, dHG, dH2G, max_iter=100, tol=1e-6):
    d_obs = np.zeros_like(H0)
    epsilon = 1e-10
    for i in range(len(H0)):
        H0_i, G0_i = H0[i], G0[i]
        G_free = G0_i
        for _ in range(max_iter):
            aH = 2 * Kd + 2 * KH2G * G_free
            bH = 1 + KHG * G_free
            cH = -H0_i
            rootsH = np.roots([aH, bH, cH])
            realH = rootsH[np.isreal(rootsH)].real
            H_free = np.min(realH[realH > epsilon]) if len(realH[realH > epsilon]) > 0 else epsilon

            aG = KH2G * H_free**2 + KHG * H_free + 1
            bG = -G0_i
            rootsG = np.roots([aG, bG])
            realG = rootsG[np.isreal(rootsG)].real
            G_free_new = np.min(realG[realG > epsilon]) if len(realG[realG > epsilon]) > 0 else epsilon

            if np.abs(G_free_new - G_free) < tol:
                break
            G_free = G_free_new

        numerator = dG * Kd * H_free**2 + dHG * KHG * H_free * G_free + dH2G * KH2G * H_free**2 * G_free
        denominator = Kd * H_free**2 + KHG * H_free * G_free + KH2G * H_free**2 * G_free
        d_obs[i] = numerator / denominator
    return d_obs

def model_definitions(H0, G0, d_delta_exp):
    return {
        "1:1": {
            "function": binding_isotherm_1_1,
            "initial_guess": [100, 100, 100],
            "bounds": ([0, -np.inf, -np.inf], [np.inf, np.inf, np.inf]),
            "lambda": lambda H0, Ka, dG, dHG: binding_isotherm_1_1(H0, G0, Ka, dG, dHG)
        },
        "1:2": {
            "function": binding_isotherm_1_2,
            "initial_guess": [100, 100, 100, 100, 100],
            "bounds": ([0, 0, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf, np.inf]),
            "lambda": lambda H0, K1, K2, dG, dHG, dHG2: binding_isotherm_1_2(H0, G0, K1, K2, dG, dHG, dHG2)
        },
        "2:1": {
            "function": binding_isotherm_2_1,
            "initial_guess": [100, 100, 100, 100, 100],
            "bounds": ([0, 0, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf, np.inf]),
            "lambda": lambda H0, K1, K2, dG, dHG, dH2G: binding_isotherm_2_1(H0, G0, K1, K2, dG, dHG, dH2G)
        },
        "dimer": {
            "function": binding_dimer,
            "initial_guess": [100, 100, 100, 100],
            "bounds": ([0, 0, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf]),
            "lambda": lambda H0, Ka, Kd, d_inf_1, d_inf_2: binding_dimer(H0, G0, Ka, Kd, d_inf_1, d_inf_2)
        },
        "multi": {
            "function": multi_model,
            "initial_guess": [100, 100, 100, 100, 100, 100],
            "bounds": ([0, 0, 0, -np.inf, -np.inf, -np.inf], [np.inf]*6),
            "lambda": lambda H0, KHG, Kd, KH2G, dG, dHG, dH2G: multi_model(H0, G0, KHG, Kd, KH2G, dG, dHG, dH2G)
        }
    }
