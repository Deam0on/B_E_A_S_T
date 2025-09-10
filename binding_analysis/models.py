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
        b = Ka * (2 * Kd * H0_i - Kd * G0_i + 1)
        c = Ka * (H0_i - G0_i) + 1
        d = -G0_i

        coefficients = [a, b, c, d]
        roots = np.roots(coefficients)
        real_roots = roots[np.isreal(roots)].real
        positive_real_roots = real_roots[real_roots > epsilon]

        G_free = (
            np.min(positive_real_roots) if len(positive_real_roots) > 0 else epsilon
        )

        # Calculate complex concentrations
        denominator = 1 + Ka * G_free + Ka * Kd * G_free**2
        HG = (Ka * H0_i * G_free) / denominator
        HG2 = (Ka * Kd * H0_i * G_free**2) / denominator

        # Calculate guest molecules in each environment
        G_in_HG = HG        # 1 guest per HG complex
        G_in_HG2 = 2 * HG2  # 2 guests per HG₂ complex
        
        # Standard weighted average: δ_obs = Σ(δ_i × [species_i]) / [total_guest]
        # For Δδ: free guest contributes 0 (reference point)
        d_delta_comp[i] = (d_inf_1 * G_in_HG + d_inf_2 * G_in_HG2) / G0_i

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
        b = Ka * (2 * Kd * G0_i - Kd * H0_i + 1)
        c = Ka * (G0_i - H0_i) + 1
        d = -H0_i

        coefficients = [a, b, c, d]
        roots = np.roots(coefficients)
        real_roots = roots[np.isreal(roots)].real
        positive_real_roots = real_roots[real_roots > epsilon]

        H_free = (
            np.min(positive_real_roots) if len(positive_real_roots) > 0 else epsilon
        )

        # Calculate complex concentrations
        denominator = 1 + Ka * H_free + Ka * Kd * H_free**2
        HG = (Ka * H_free * G0_i) / denominator
        H2G = (Ka * Kd * H_free**2 * G0_i) / denominator

        # For 2:1 binding, each guest is in exactly one environment
        # HG contains [HG] guest molecules
        # H2G contains [H2G] guest molecules
        # Standard weighted average: δ_obs = Σ(δ_i × [species_i]) / [total_guest]
        d_delta_comp[i] = (d_inf_1 * HG + d_inf_2 * H2G) / G0_i

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
        d_inf_1: Chemical shift change for guest in HG complex (guest-observed)
        d_inf_2: Chemical shift change for free guest (usually 0, for reference)

    Returns:
        Calculated chemical shift changes (Δδ) for guest-observed NMR
    """
    H0, G0 = np.asarray(H0), np.asarray(G0)
    d_delta_comp = np.zeros_like(H0, dtype=float)
    epsilon = 1e-10

    for i in range(len(H0)):
        H0_i, G0_i = H0[i], G0[i]

        # Cubic equation coefficients for H_free
        # Derived from: 2*Kd*Ka*[H]³ + 2*Kd*[H]² + Ka*(G0-H0)*[H] + H0 = 0
        a = 2 * Kd * Ka
        b = 2 * Kd  
        c = Ka * (G0_i - H0_i)
        d = H0_i

        # Handle degenerate cases
        if abs(a) < epsilon:  # No cubic term
            if abs(b) < epsilon:  # Linear equation
                if abs(c) < epsilon:
                    H_free = epsilon
                else:
                    H_free = -d / c
            else:  # Quadratic equation
                discriminant = c**2 - 4*b*d
                if discriminant >= 0:
                    sqrt_disc = np.sqrt(discriminant)
                    root1 = (-c + sqrt_disc) / (2*b)
                    root2 = (-c - sqrt_disc) / (2*b)
                    H_free = max(root1, root2) if max(root1, root2) > epsilon else epsilon
                else:
                    H_free = epsilon
        else:  # Full cubic equation
            coefficients = [a, b, c, d]
            roots = np.roots(coefficients)
            real_roots = roots[np.isreal(roots)].real
            positive_real_roots = real_roots[real_roots > epsilon]
            
            if len(positive_real_roots) > 0:
                # Choose the physically meaningful root
                # Usually the smallest positive root
                H_free = np.min(positive_real_roots)
            else:
                H_free = epsilon

        # Calculate species concentrations
        G_free = G0_i / (1 + Ka * H_free)
        HG = Ka * H_free * G_free
        H2 = Kd * H_free**2

        # For guest-observed NMR: weighted average of guest environments
        # G_free contributes 0 (reference), HG contributes d_inf_1
        d_delta_comp[i] = (d_inf_1 * HG) / G0_i

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
        dG: Chemical shift of free guest (reference, usually 0)
        dHG: Chemical shift of guest in HG complex
        dH2G: Chemical shift of guest in H₂G complex
        max_iter: Maximum iterations for convergence
        tol: Convergence tolerance

    Returns:
        Calculated chemical shift changes (Δδ) for guest-observed NMR
    """
    H0, G0 = np.asarray(H0), np.asarray(G0)
    d_delta_comp = np.zeros_like(H0, dtype=float)
    epsilon = 1e-10

    for i in range(len(H0)):
        H0_i, G0_i = H0[i], G0[i]
        
        # Initial guess for free concentrations
        H_free = H0_i / 2  # Start with reasonable guess
        G_free = G0_i / 2
        
        # Iterative solution for equilibrium concentrations
        for iteration in range(max_iter):
            H_free_old = H_free
            G_free_old = G_free
            
            # Calculate complex concentrations from current free concentrations
            H2 = Kd * H_free**2
            HG = KHG * H_free * G_free
            H2G = KH2G * H2 * G_free
            
            # Update free concentrations from mass balance
            # From host mass balance: H0 = H_free + HG + 2*H2 + 2*H2G
            # From guest mass balance: G0 = G_free + HG + H2G
            
            # Guest mass balance (simpler)
            G_free = G0_i / (1 + KHG * H_free + KH2G * Kd * H_free**2)
            
            # Host mass balance (quadratic in H_free)
            # H0 = H_free + KHG*H_free*G_free + 2*Kd*H_free^2 + 2*KH2G*Kd*H_free^2*G_free
            # H0 = H_free * (1 + KHG*G_free + 2*Kd*H_free + 2*KH2G*Kd*H_free*G_free)
            
            # Rearrange to: 2*Kd*(1 + KH2G*G_free)*H_free^2 + (1 + KHG*G_free)*H_free - H0 = 0
            a = 2 * Kd * (1 + KH2G * G_free)
            b = 1 + KHG * G_free
            c = -H0_i
            
            if abs(a) < epsilon:  # Linear case
                H_free = -c / b if abs(b) > epsilon else epsilon
            else:  # Quadratic case
                discriminant = b**2 - 4*a*c
                if discriminant >= 0:
                    sqrt_disc = np.sqrt(discriminant)
                    root1 = (-b + sqrt_disc) / (2*a)
                    root2 = (-b - sqrt_disc) / (2*a)
                    H_free = max(root1, root2) if max(root1, root2) > epsilon else epsilon
                else:
                    H_free = epsilon
            
            # Check convergence
            if (abs(H_free - H_free_old) < tol and abs(G_free - G_free_old) < tol):
                break
        
        # Calculate final species concentrations
        H2 = Kd * H_free**2
        HG = KHG * H_free * G_free
        H2G = KH2G * H2 * G_free
        
        # Guest-observed chemical shift calculation
        # Only species containing guest contribute: G_free, HG, H2G
        # Standard weighted average: δ_obs = Σ(δ_i × [species_i]) / [total_guest]
        
        total_guest = G_free + HG + H2G
        if total_guest > epsilon:
            delta_obs = (dG * G_free + dHG * HG + dH2G * H2G) / total_guest
        else:
            delta_obs = dG  # Fallback to free guest shift
        
        # Convert to Δδ relative to free guest (dG as reference)
        d_delta_comp[i] = delta_obs - dG

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
