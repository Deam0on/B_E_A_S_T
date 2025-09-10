# Mathematical Derivations of Binding Isotherm Models

This document provides detailed mathematical derivations for the 1:1, 1:2, and 2:1 host-guest binding models used in BEAST (Binding Evaluation and Analysis Software Tool).

## Table of Contents
- [1:1 Host-Guest Binding Model](#11-host-guest-binding-model)
- [1:2 Host-Guest Binding Model](#12-host-guest-binding-model)
- [2:1 Host-Guest Binding Model](#21-host-guest-binding-model)
- [Chemical Shift Calculations](#chemical-shift-calculations)
- [References](#references)

---

## 1:1 Host-Guest Binding Model

### Equilibrium Expression

The simplest binding model describes the formation of a 1:1 complex between host (H) and guest (G):

```
H + G ⇌ HG
```

**Equilibrium constant:**

$$K_a = \frac{[HG]}{[H][G]}$$

### Mass Balance Equations

**Total host concentration:**

$$[H_0] = [H] + [HG]$$

**Total guest concentration:**

$$[G_0] = [G] + [HG]$$

### Algebraic Solution

From the mass balance equations:

$$[H] = [H_0] - [HG]$$
$$[G] = [G_0] - [HG]$$

Substituting into the equilibrium expression:

$$K_a = \frac{[HG]}{([H_0] - [HG])([G_0] - [HG])}$$

Rearranging:

$$K_a([H_0] - [HG])([G_0] - [HG]) = [HG]$$

Expanding:

$$K_a([H_0][G_0] - [H_0][HG] - [G_0][HG] + [HG]^2) = [HG]$$

$$K_a[H_0][G_0] - K_a[H_0][HG] - K_a[G_0][HG] + K_a[HG]^2 = [HG]$$

Rearranging into standard quadratic form:

$$K_a[HG]^2 - (K_a[H_0] + K_a[G_0] + 1)[HG] + K_a[H_0][G_0] = 0$$

Dividing by $K_a$:

$$[HG]^2 - ([H_0] + [G_0] + \frac{1}{K_a})[HG] + [H_0][G_0] = 0$$

**Quadratic formula solution:**

$$[HG] = \frac{([H_0] + [G_0] + \frac{1}{K_a}) \pm \sqrt{([H_0] + [G_0] + \frac{1}{K_a})^2 - 4[H_0][G_0]}}{2}$$

---

## 1:2 Host-Guest Binding Model

### Equilibrium Expressions

This model describes sequential binding of two guests to one host:

```
H + G ⇌ HG     (Ka)
HG + G ⇌ HG₂   (Kd)
```

**Equilibrium constants:**

$$K_a = \frac{[HG]}{[H][G]}$$
$$K_d = \frac{[HG_2]}{[HG][G]}$$

### Mass Balance Equations

**Total host concentration:**

$$[H_0] = [H] + [HG] + [HG_2]$$

**Total guest concentration:**

$$[G_0] = [G] + [HG] + 2[HG_2]$$

### Algebraic Solution

From the equilibrium expressions:

$$[HG] = K_a[H][G]$$
$$[HG_2] = K_d[HG][G] = K_a K_d[H][G]^2$$

From host mass balance:

$$[H_0] = [H] + K_a[H][G] + K_a K_d[H][G]^2$$
$$[H_0] = [H](1 + K_a[G] + K_a K_d[G]^2)$$

Therefore:

$$[H] = \frac{[H_0]}{1 + K_a[G] + K_a K_d[G]^2}$$

Substituting back:

$$[HG] = \frac{K_a[H_0][G]}{1 + K_a[G] + K_a K_d[G]^2}$$
$$[HG_2] = \frac{K_a K_d[H_0][G]^2}{1 + K_a[G] + K_a K_d[G]^2}$$

From guest mass balance:

$$[G_0] = [G] + \frac{K_a[H_0][G]}{1 + K_a[G] + K_a K_d[G]^2} + 2\frac{K_a K_d[H_0][G]^2}{1 + K_a[G] + K_a K_d[G]^2}$$

Multiplying through by the denominator:

$$[G_0](1 + K_a[G] + K_a K_d[G]^2) = [G](1 + K_a[G] + K_a K_d[G]^2) + K_a[H_0][G] + 2K_a K_d[H_0][G]^2$$

Rearranging:

$$[G_0] + [G_0]K_a[G] + [G_0]K_a K_d[G]^2 = [G] + K_a[G]^2 + K_a K_d[G]^3 + K_a[H_0][G] + 2K_a K_d[H_0][G]^2$$

Collecting terms:

$$K_a K_d[G]^3 + K_a(2K_d[H_0] - K_d[G_0] + 1)[G]^2 + K_a([H_0] - [G_0])[G] + ([G_0]) = 0$$

**Cubic equation in [G] (free guest concentration):**

$$a[G]^3 + b[G]^2 + c[G] + d = 0$$

Where:

- $a = K_a K_d$
- $b = K_a(2K_d[H_0] - K_d[G_0] + 1)$
- $c = K_a([H_0] - [G_0]) + 1$
- $d = -[G_0]$

---

## 2:1 Host-Guest Binding Model

### Equilibrium Expressions

This model describes sequential binding of two hosts to one guest:

```
H + G ⇌ HG      (Ka)
H + HG ⇌ H₂G    (Kd)
```

**Equilibrium constants:**

$$K_a = \frac{[HG]}{[H][G]}$$
$$K_d = \frac{[H_2G]}{[H][HG]}$$

### Mass Balance Equations

**Total host concentration:**

$$[H_0] = [H] + [HG] + 2[H_2G]$$

**Total guest concentration:**

$$[G_0] = [G] + [HG] + [H_2G]$$

### Algebraic Solution

From the equilibrium expressions:

$$[HG] = K_a[H][G]$$
$$[H_2G] = K_d[H][HG] = K_a K_d[H]^2[G]$$

From guest mass balance:

$$[G_0] = [G] + K_a[H][G] + K_a K_d[H]^2[G]$$
$$[G_0] = [G](1 + K_a[H] + K_a K_d[H]^2)$$

Therefore:

$$[G] = \frac{[G_0]}{1 + K_a[H] + K_a K_d[H]^2}$$

Substituting back:

$$[HG] = \frac{K_a[H][G_0]}{1 + K_a[H] + K_a K_d[H]^2}$$
$$[H_2G] = \frac{K_a K_d[H]^2[G_0]}{1 + K_a[H] + K_a K_d[H]^2}$$

From host mass balance:

$$[H_0] = [H] + \frac{K_a[H][G_0]}{1 + K_a[H] + K_a K_d[H]^2} + 2\frac{K_a K_d[H]^2[G_0]}{1 + K_a[H] + K_a K_d[H]^2}$$

Multiplying through by the denominator:

$$[H_0](1 + K_a[H] + K_a K_d[H]^2) = [H](1 + K_a[H] + K_a K_d[H]^2) + K_a[H][G_0] + 2K_a K_d[H]^2[G_0]$$

Rearranging:

$$[H_0] + [H_0]K_a[H] + [H_0]K_a K_d[H]^2 = [H] + K_a[H]^2 + K_a K_d[H]^3 + K_a[G_0][H] + 2K_a K_d[G_0][H]^2$$

Collecting terms:

$$K_a K_d[H]^3 + K_a(2K_d[G_0] - K_d[H_0] + 1)[H]^2 + K_a([G_0] - [H_0])[H] + ([H_0]) = 0$$

**Cubic equation in [H] (free host concentration):**

$$a[H]^3 + b[H]^2 + c[H] + d = 0$$

Where:

- $a = K_a K_d$
- $b = K_a(2K_d[G_0] - K_d[H_0] + 1)$
- $c = K_a([G_0] - [H_0]) + 1$
- $d = -[H_0]$

---

## Chemical Shift Calculations

### Fundamental Principle

The observed chemical shift in NMR is the population-weighted average of all species:

$$\delta_{obs} = \frac{\sum_i \delta_i [species_i]}{[total]}$$

For guest-observed experiments:

$$\delta_{obs} = \frac{\delta_{free}[G_{free}] + \delta_{HG}[G_{in\_HG}] + \delta_{complex}[G_{in\_complex}]}{[G_0]}$$

### Delta-Delta (Δδ) Calculation

Chemical shift changes are calculated relative to the free state:

$$\Delta\delta = \delta_{obs} - \delta_{free}$$

Since $\delta_{free}$ serves as the reference (zero point):

$$\Delta\delta = \frac{\delta_{HG}[G_{in\_HG}] + \delta_{complex}[G_{in\_complex}]}{[G_0]}$$

### Model-Specific Implementations

#### 1:1 Model
- Guest molecules: $[G_{\text{free}}]$ and $[G_{\text{in HG}}] = [HG]$
- $\Delta\delta = \frac{\delta_{\text{inf}} \cdot [HG]}{[G_0]}$

#### 1:2 Model
- Guest molecules: $[G_{\text{free}}]$, $[G_{\text{in HG}}] = [HG]$, and $[G_{\text{in HG}_2}] = 2[HG_2]$
- $\Delta\delta = \frac{\delta_{\text{inf1}} \cdot [HG] + \delta_{\text{inf2}} \cdot 2[HG_2]}{[G_0]}$

#### 2:1 Model
- Guest molecules: $[G_{\text{free}}]$, $[G_{\text{in HG}}] = [HG]$, and $[G_{\text{in H}_2\text{G}}] = [H_2G]$
- $\Delta\delta = \frac{\delta_{\text{inf1}} \cdot [HG] + \delta_{\text{inf2}} \cdot [H_2G]}{[G_0]}$

### Parameter Interpretation

- **$\delta_{inf}$, $\delta_{inf1}$, $\delta_{inf2}$**: These represent the chemical shift changes of the guest when fully bound in the respective complex environments
- **Physical meaning**: These are the limiting chemical shifts that would be observed if all guest molecules were in that specific binding environment

---

## References

1. **Fielding, L.** *Determination of association constants (Ka) from solution NMR data.* Tetrahedron 2000, 56, 6151-6170.

2. **Hirose, K.** *A practical guide for the determination of binding constants.* J. Incl. Phenom. Macrocycl. Chem. 2001, 39, 193-209.

3. **Connors, K.A.** *Binding Constants: The Measurement of Molecular Complex Stability.* John Wiley & Sons: New York, 1987.

4. **Thordarson, P.** *Determining association constants from titration experiments in supramolecular chemistry.* Chem. Soc. Rev. 2011, 40, 1305-1323.

5. **Ulatowski, F.; Dąbrowa, K.; Bałakier, T.; Jurczak, J.** *Recognizing the Limited Applicability of Job Plots in Studying Host–Guest Interactions in Supramolecular Chemistry.* J. Org. Chem. 2016, 81, 1746-1756.

---

*This document provides the mathematical foundation for the binding isotherm models implemented in BEAST. The derivations follow standard approaches in supramolecular chemistry and NMR titration analysis.*