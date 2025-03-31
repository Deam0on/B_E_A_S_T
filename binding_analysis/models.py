
import numpy as np

def binding_isotherm_1_1(H0, G0, d_delta_exp, Ka, d_inf):
    term = G0 + H0 + (1 / Ka)
    sqrt_term = np.sqrt(term * term - 4 * G0 * H0)
    HG = 0.5 * (term - sqrt_term)
    d_free = d_delta_exp[0]
    d_delta_comp = (d_free * (HG - G0) / G0) + (d_inf * (HG / G0))
    return d_delta_comp
    
def binding_isotherm_1_2(H0, G0, Ka, Kd, d_inf_1, d_inf_2):
    d_delta_comp = np.zeros_like(H0)
    for i in range(len(H0)):
        epsilon = 1e-10
        H0_i = H0[i]
        G0_i = G0[i]
        a = (Ka * Kd)
        b = (Ka * ((2 * Kd * H0_i) - (Kd * G0_i) + 1))
        c = (Ka * (H0_i - G0_i) + 1)
        d = -1 * G0_i
        coefficients = [a, b, c, d]
        roots = np.roots(coefficients)
        real_roots = roots[np.isreal(roots)].real
        positive_real_roots = real_roots[real_roots > epsilon]
        G_free = np.min(positive_real_roots) if len(positive_real_roots) > 0 else epsilon
        d_delta_comp[i] = ((d_inf_1 * G0_i * Ka * G_free) + (d_inf_2 * G0_i * Ka * Kd * G_free**2)) / (1 + (Ka * G_free) + (Ka * Kd * G_free**2))
    return d_delta_comp

def binding_isotherm_2_1(H0, G0, Ka, Kd, d_inf_1, d_inf_2):
    d_delta_comp = np.zeros_like(H0)
    for i in range(len(H0)):
        epsilon = 1e-10
        H0_i = H0[i]
        G0_i = G0[i]
        a = (Ka * Kd)
        b = (Ka * ((2 * Kd * G0_i) - (Kd * H0_i) + 1))
        c = (Ka * (G0_i - H0_i) + 1)
        d = -H0_i
        coefficients = [a, b, c, d]
        roots = np.roots(coefficients)
        real_roots = roots[np.isreal(roots)].real
        positive_real_roots = real_roots[real_roots > epsilon]
        H_free = np.min(positive_real_roots) if len(positive_real_roots) > 0 else epsilon
        d_delta_comp[i] = ((d_inf_1 * G0_i * Ka * H_free) + (2 * d_inf_2 * G0_i * Ka * Kd * H_free**2)) / (G0_i * (1 + (Ka * H_free) + (Ka * Kd* H_free**2)))
    return d_delta_comp


def binding_dimer(H0, G0, Ka, Kd, d_inf_1, d_inf_2):
    d_obs = np.zeros_like(H0)
    epsilon = 1e-10

    for i in range(len(H0)):
        H0_i = H0[i]
        G0_i = G0[i]

        a = -2 * Kd
        b = -1 * (G0_i * Ka + 1)
        c = H0_i
        roots = np.roots([a, b, c])
        real_roots = roots[np.isreal(roots)].real
        H_free = np.min(real_roots[real_roots > epsilon]) if len(real_roots[real_roots > epsilon]) > 0 else epsilon

        numerator = (d_inf_1 * G0_i) / (1 + Ka * H_free) + (d_inf_2 * Ka * H_free * G0_i)
        denominator = (G0_i / (1 + Ka * H_free)) + (Ka * H_free * G0_i)
        d_obs[i] = numerator / denominator

    # Convert to Δδ
    d_delta_comp = d_obs - d_obs[0]
    d_delta_comp[0] = 0
    return d_delta_comp

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

    # Convert to Δδ
    d_delta_comp = d_obs - d_obs[0]
    d_delta_comp[0] = 0
    return d_delta_comp

def model_definitions(H0, G0, d_delta_exp):
    return {
        "1:1": {
            "function": binding_isotherm_1_1,
            "initial_guess": [100, 100, 100],  # Ka, d_free, d_inf
            "bounds": ([0, -np.inf, -np.inf], [np.inf, np.inf, np.inf]),
            "lambda": lambda H0, Ka, d_free, d_inf: binding_isotherm_1_1(H0, G0, d_delta_exp, Ka, d_inf)
        },
        "1:2": {
            "function": binding_isotherm_1_2,
            "initial_guess": [100, 100, 100, 100],
            "bounds": ([0, 0, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf]),
            "lambda": lambda H0, Ka, Kd, d_inf_1, d_inf_2: binding_isotherm_1_2(H0, G0, Ka, Kd, d_inf_1, d_inf_2)
        },
        "2:1": {
            "function": binding_isotherm_2_1,
            "initial_guess": [100, 100, 100, 100],
            "bounds": ([0, 0, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf]),
            "lambda": lambda H0, Ka, Kd, d_inf_1, d_inf_2: binding_isotherm_2_1(H0, G0, Ka, Kd, d_inf_1, d_inf_2)
        },
        "dimer": {
            "function": binding_dimer,
            "initial_guess": [100, 100, 100, 100],
            "bounds": ([0, 0, -np.inf, -np.inf], [np.inf]*4),
            "lambda": lambda H0, Ka, Kd, d_inf_1, d_inf_2: binding_dimer(H0, G0, Ka, Kd, d_inf_1, d_inf_2)
        },
        "multi": {
            "function": multi_model,
            "initial_guess": [100, 100, 100, 100, 100, 100],
            "bounds": ([0, 0, 0, -np.inf, -np.inf, -np.inf], [np.inf]*6),
            "lambda": lambda H0, KHG, Kd, KH2G, dG, dHG, dH2G: multi_model(H0, G0, KHG, Kd, KH2G, dG, dHG, dH2G)
        }
    }
