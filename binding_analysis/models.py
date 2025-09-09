"""
Binding isotherm models for NMR titration data analysis.

This module contains mathematical models for various host-guest binding scenarios
including 1:1, 1:2, 2:1 binding, dimerization, and multi-equilibrium systems.

Author: Filip Hládek
License: MIT
"""

from typing import Any, Callable, Dict, Tuple, Union

import numpy as np

# Type aliases for clarity
ArrayLike = Union[np.ndarray, list, float]


def binding_isotherm_1_1(
    H0: ArrayLike, G0: ArrayLike, Ka: float, d_free: float, d_inf: float
) -> np.ndarray:
    """
    Calculate chemical shift changes for 1:1 host-guest binding.

    This model describes simple 1:1 binding between host (H) and guest (G):
    H + G ⇌ HG

    Args:
        H0: Total host concentration(s)
        G0: Total guest concentration(s)
        Ka: Association constant (M⁻¹)
        d_free: Chemical shift of free host
        d_inf: Chemical shift of bound host (HG complex)

    Returns:
        Calculated chemical shift changes (Δδ)
    """
    H0, G0 = np.asarray(H0), np.asarray(G0)

    # Quadratic formula solution for HG concentration
    term = G0 + H0 + (1 / Ka)
    sqrt_term = np.sqrt(term * term - 4 * G0 * H0)
    HG = 0.5 * (term - sqrt_term)

    # Chemical shift calculation
    # d_delta_comp = (d_free * (G0-HG) / G0) + (d_inf * (HG / G0))
    # fixed to be delta-delta
    d_delta_comp = (d_inf * (HG / G0))
    return d_delta_comp


def binding_isotherm_1_2(
    H0: ArrayLike, G0: ArrayLike, Ka: float, Kd: float, d_inf_1: float, d_inf_2: float
) -> np.ndarray:
    """
    Calculate chemical shift changes for 1:2 host-guest binding.

    This model describes sequential binding of two guests to one host:
    H + G ⇌ HG    (Ka)
    HG + G ⇌ HG₂  (Kd)

    Args:
        H0: Total host concentration(s)
        G0: Total guest concentration(s)
        Ka: First binding constant (M⁻¹)
        Kd: Second binding constant (M⁻¹)
        d_inf_1: Chemical shift contribution from HG
        d_inf_2: Chemical shift contribution from HG₂

    Returns:
        Calculated chemical shift changes (Δδ)
    """
    H0, G0 = np.asarray(H0), np.asarray(G0)
    d_delta_comp = np.zeros_like(H0, dtype=float)
    epsilon = 1e-10

    for i in range(len(H0)):
        H0_i, G0_i = H0[i], G0[i]

        # Cubic equation coefficients for G_free
        a = Ka * Kd
        b = Ka * ((2 * Kd * H0_i) - (Kd * G0_i) + 1)
        c = Ka * (H0_i - G0_i) + 1
        d = -G0_i

        coefficients = [a, b, c, d]
        roots = np.roots(coefficients)
        real_roots = roots[np.isreal(roots)].real
        positive_real_roots = real_roots[real_roots > epsilon]

        G_free = (
            np.min(positive_real_roots) if len(positive_real_roots) > 0 else epsilon
        )

        # Chemical shift calculation
        numerator = d_inf_1 * G0_i * Ka * G_free + d_inf_2 * G0_i * Ka * Kd * G_free**2
        denominator = 1 + Ka * G_free + Ka * Kd * G_free**2

        d_delta_comp[i] = numerator / denominator

    return d_delta_comp


def binding_isotherm_2_1(
    H0: ArrayLike, G0: ArrayLike, Ka: float, Kd: float, d_inf_1: float, d_inf_2: float
) -> np.ndarray:
    """
    Calculate chemical shift changes for 2:1 host-guest binding.

    This model describes sequential binding of two hosts to one guest:
    H + G ⇌ HG     (Ka)
    H + HG ⇌ H₂G   (Kd)

    Args:
        H0: Total host concentration(s)
        G0: Total guest concentration(s)
        Ka: First binding constant (M⁻¹)
        Kd: Second binding constant (M⁻¹)
        d_inf_1: Chemical shift contribution from HG
        d_inf_2: Chemical shift contribution from H₂G

    Returns:
        Calculated chemical shift changes (Δδ)
    """
    H0, G0 = np.asarray(H0), np.asarray(G0)
    d_delta_comp = np.zeros_like(H0, dtype=float)
    epsilon = 1e-10

    for i in range(len(H0)):
        H0_i, G0_i = H0[i], G0[i]

        # Cubic equation coefficients for H_free
        a = Ka * Kd
        b = Ka * ((2 * Kd * G0_i) - (Kd * H0_i) + 1)
        c = Ka * (G0_i - H0_i) + 1
        d = -H0_i

        coefficients = [a, b, c, d]
        roots = np.roots(coefficients)
        real_roots = roots[np.isreal(roots)].real
        positive_real_roots = real_roots[real_roots > epsilon]

        H_free = (
            np.min(positive_real_roots) if len(positive_real_roots) > 0 else epsilon
        )

        # Chemical shift calculation
        numerator = (
            d_inf_1 * H0_i * Ka * H_free + 2 * d_inf_2 * H0_i * Ka * Kd * H_free**2
        )
        if G0_i <= 0:
            d_delta_comp[i] = 0
        else:
            denominator = G0_i * (1 + Ka * H_free + Ka * Kd * H_free**2)
            d_delta_comp[i] = numerator / denominator

        # d_delta_comp[i] = numerator / denominator

    return d_delta_comp


def binding_dimer(
    H0: ArrayLike, G0: ArrayLike, Ka: float, Kd: float, d_inf_1: float, d_inf_2: float
) -> np.ndarray:
    """
    Calculate chemical shift changes for host-guest binding with host dimerization.

    This model describes binding with simultaneous host dimerization:
    H + G ⇌ HG     (Ka)
    2H ⇌ H₂       (Kd)

    Args:
        H0: Total host concentration(s)
        G0: Total guest concentration(s)
        Ka: Host-guest binding constant (M⁻¹)
        Kd: Host dimerization constant (M⁻¹)
        d_inf_1: Chemical shift of free host
        d_inf_2: Chemical shift of bound host

    Returns:
        Calculated chemical shift changes (Δδ) relative to first point
    """
    H0, G0 = np.asarray(H0), np.asarray(G0)
    d_obs = np.zeros_like(H0, dtype=float)
    epsilon = 1e-10

    for i in range(len(H0)):
        H0_i, G0_i = H0[i], G0[i]

        # Quadratic equation for H_free
        a = -2 * Kd
        b = -(G0_i * Ka + 1)
        c = H0_i

        roots = np.roots([a, b, c])
        real_roots = roots[np.isreal(roots)].real
        valid_roots = real_roots[real_roots > epsilon]

        H_free = np.min(valid_roots) if len(valid_roots) > 0 else epsilon

        # Chemical shift calculation
        numerator = (d_inf_1 * G0_i) / (1 + Ka * H_free) + (
            d_inf_2 * Ka * H_free * G0_i / (1 + Ka * H_free)
        )
        denominator = (G0_i / (1 + Ka * H_free)) + (Ka * H_free * G0_i / (1 + Ka * H_free))

        d_obs[i] = numerator / denominator

    # Convert to delta-delta (relative to first point)
    d_delta_comp = d_obs - d_obs[0]
    d_delta_comp[0] = 0
    return d_delta_comp


def multi_model(
    H0: ArrayLike,
    G0: ArrayLike,
    KHG: float,
    Kd: float,
    KH2G: float,
    dG: float,
    dHG: float,
    dH2G: float,
    max_iter: int = 100000,
    tol: float = 1e-6,
) -> np.ndarray:
    """
    Calculate chemical shift changes for multi-equilibrium system.

    This model describes a complex system with multiple equilibria:
    H + G ⇌ HG      (KHG)
    2H ⇌ H₂        (Kd)
    H₂ + G ⇌ H₂G    (KH2G)

    Args:
        H0: Total host concentration(s)
        G0: Total guest concentration(s)
        KHG: Host-guest binding constant (M⁻¹)
        Kd: Host dimerization constant (M⁻¹)
        KH2G: Dimer-guest binding constant (M⁻¹)
        dG: Chemical shift of free guest
        dHG: Chemical shift of HG complex
        dH2G: Chemical shift of H₂G complex
        max_iter: Maximum iterations for convergence
        tol: Convergence tolerance

    Returns:
        Calculated chemical shift changes (Δδ) relative to first point
    """
    H0, G0 = np.asarray(H0), np.asarray(G0)
    d_obs = np.zeros_like(H0, dtype=float)
    epsilon = 1e-10

    for i in range(len(H0)):
        H0_i, G0_i = H0[i], G0[i]
        G_free = G0_i  # Initial guess

        # Iterative solution for equilibrium concentrations
        for iteration in range(max_iter):
            # Solve for H_free given G_free
            aH = 2 * Kd + 2 * KH2G * G_free
            bH = 1 + KHG * G_free
            cH = -H0_i

            rootsH = np.roots([aH, bH, cH])
            realH = rootsH[np.isreal(rootsH)].real
            valid_H = realH[realH > epsilon]

            H_free = np.min(valid_H) if len(valid_H) > 0 else epsilon

            # Solve for G_free given H_free
            aG = KH2G * H_free**2 + KHG * H_free + 1
            bG = -G0_i

            rootsG = np.roots([aG, bG])
            realG = rootsG[np.isreal(rootsG)].real
            valid_G = realG[realG > epsilon]

            G_free_new = np.min(valid_G) if len(valid_G) > 0 else epsilon

            # Check convergence
            if abs(G_free_new - G_free) < tol:
                break
            G_free = G_free_new

        # Calculate chemical shift
        numerator = (
            dG * Kd * H_free**2
            + dHG * KHG * H_free * G_free
            + dH2G * KH2G * H_free**2 * G_free
        )
        denominator = Kd * H_free**2 + KHG * H_free * G_free + KH2G * H_free**2 * G_free

        d_obs[i] = numerator / denominator if denominator > epsilon else 0

    # Convert to delta-delta (relative to first point)
    d_delta_comp = d_obs - d_obs[0]
    d_delta_comp[0] = 0
    return d_delta_comp


def model_definitions(
    H0: ArrayLike, G0: ArrayLike, d_delta_exp: ArrayLike
) -> Dict[str, Dict[str, Any]]:
    """
    Define available binding models with their parameters and constraints.

    Args:
        H0: Total host concentration(s)
        G0: Total guest concentration(s)
        d_delta_exp: Experimental chemical shift changes (not used in definitions)

    Returns:
        Dictionary of model definitions with functions, initial guesses, bounds, and lambdas
    """
    return {
        "HG": {
            "function": binding_isotherm_1_1,
            "initial_guess": [100, 100, 100],  # Ka, d_free, d_inf
            "bounds": ([0, -np.inf, -np.inf], [np.inf, np.inf, np.inf]),
            "lambda": lambda H0, Ka, d_free, d_inf: binding_isotherm_1_1(
                H0, G0, Ka, d_free, d_inf
            ),
            "description": "1:1 Host-Guest binding (H + G ⇌ HG)",
        },
        "HG₂": {
            "function": binding_isotherm_1_2,
            "initial_guess": [100, 100, 100, 100],  # Ka, Kd, d_inf_1, d_inf_2
            "bounds": ([0, 0, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf]),
            "lambda": lambda H0, Ka, Kd, d_inf_1, d_inf_2: binding_isotherm_1_2(
                H0, G0, Ka, Kd, d_inf_1, d_inf_2
            ),
            "description": "1:2 Host-Guest binding (H + G ⇌ HG, HG + G ⇌ HG₂)",
        },
        "H₂G": {
            "function": binding_isotherm_2_1,
            "initial_guess": [100, 100, 100, 100],  # Ka, Kd, d_inf_1, d_inf_2
            "bounds": ([0, 0, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf]),
            "lambda": lambda H0, Ka, Kd, d_inf_1, d_inf_2: binding_isotherm_2_1(
                H0, G0, Ka, Kd, d_inf_1, d_inf_2
            ),
            "description": "2:1 Host-Guest binding (H + G ⇌ HG, H + HG ⇌ H₂G)",
        },
        "HG + H₂": {
            "function": binding_dimer,
            "initial_guess": [100, 100, 100, 100],  # Ka, Kd, d_inf_1, d_inf_2
            "bounds": ([0, 0, -np.inf, -np.inf], [np.inf] * 4),
            "lambda": lambda H0, Ka, Kd, d_inf_1, d_inf_2: binding_dimer(
                H0, G0, Ka, Kd, d_inf_1, d_inf_2
            ),
            "description": "Host-Guest binding with dimerization (H + G ⇌ HG, 2H ⇌ H₂)",
        },
        "H₂G + HG + H₂": {
            "function": multi_model,
            "initial_guess": [
                100,
                100,
                100,
                100,
                100,
                100,
            ],  # KHG, Kd, KH2G, dG, dHG, dH2G
            "bounds": ([0, 0, 0, -np.inf, -np.inf, -np.inf], [np.inf] * 6),
            "lambda": lambda H0, KHG, Kd, KH2G, dG, dHG, dH2G: multi_model(
                H0, G0, KHG, Kd, KH2G, dG, dHG, dH2G
            ),
            "description": "Multi-equilibrium system (H + G ⇌ HG, 2H ⇌ H₂, H₂ + G ⇌ H₂G)",
        },
    }
