# Mathematical Derivations of Binding Isotherm Models

This document provides detailed mathematical derivations for the 1:1, 1:2, and 2:1 host-guest binding models used in BEAST (Binding Evaluation and Analysis Software Tool).

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

$$K_a K_d[G]^3 + K_a(2K_d[H_0] - K_d[G_0] + 1)[G]^2 + (K_a([H_0] - [G_0])+1)[G] - ([G_0]) = 0$$

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

$$K_a K_d[H]^3 + K_a(2K_d[G_0] - K_d[H_0] + 1)[H]^2 + (K_a([G_0] - [H_0])+1)[H] - ([H_0]) = 0$$

**Cubic equation in [H] (free host concentration):**

$$a[H]^3 + b[H]^2 + c[H] + d = 0$$

Where:

- $a = K_a K_d$
- $b = K_a(2K_d[G_0] - K_d[H_0] + 1)$
- $c = K_a([G_0] - [H_0]) + 1$
- $d = -[H_0]$

---

## Host Dimerization Model (HG + H₂)

### Equilibrium Expressions

This model describes host-guest binding in competition with host dimerization:

```
H + G ⇌ HG     (Ka)
2H ⇌ H₂       (Kd)
```

**Equilibrium constants:**

$$K_a = \frac{[HG]}{[H][G]}$$
$$K_d = \frac{[H_2]}{[H]^2}$$

### Mass Balance Equations

**Total host concentration:**

$$[H_0] = [H] + [HG] + 2[H_2]$$

**Total guest concentration:**

$$[G_0] = [G] + [HG]$$

### Algebraic Solution

From the equilibrium expressions:

$$[HG] = K_a[H][G]$$
$$[H_2] = K_d[H]^2$$

From guest mass balance:

$$[G_0] = [G] + K_a[H][G]$$
$$[G_0] = [G](1 + K_a[H])$$

Therefore:

$$[G] = \frac{[G_0]}{1 + K_a[H]}$$

Substituting back:

$$[HG] = \frac{K_a[H][G_0]}{1 + K_a[H]}$$

From host mass balance:

$$[H_0] = [H] + \frac{K_a[H][G_0]}{1 + K_a[H]} + 2K_d[H]^2$$

Multiplying through by $(1 + K_a[H])$:

$$[H_0](1 + K_a[H]) = [H](1 + K_a[H]) + K_a[H][G_0] + 2K_d[H]^2(1 + K_a[H])$$

Expanding:

$$[H_0] + [H_0]K_a[H] = [H] + K_a[H]^2 + K_a[H][G_0] + 2K_d[H]^2 + 2K_dK_a[H]^3$$

Rearranging:

$$2K_dK_a[H]^3 + (K_a + 2K_d)[H]^2 + (K_a([G_0] - [H_0])+1)[H] - [H_0] = 0$$

**Cubic equation in [H] (free host concentration):**

$$a[H]^3 + b[H]^2 + c[H] + d = 0$$

Where:

- $a = 2K_aK_d$
- $b = K_a + 2K_d$
- $c = K_a([G_0] - [H_0]) + 1$
- $d = -[H_0]$

---

## Multi-Equilibrium Model (HG + H₂ + H₂G)

### Equilibrium Expressions

This model describes a complex system with three competing equilibria:

```
H + G ⇌ HG      (KHG)
2H ⇌ H₂        (Kd)
H + HG ⇌ H₂G    (KH2G)
```

**Equilibrium constants:**
$$K_{HG} = \frac{[HG]}{[H][G]}$$
$$K_d = \frac{[H_2]}{[H]^2}$$
$$K_{H_2G} = \frac{[H_2G]}{[H][HG]}$$

### Mass Balance Equations

**Total host concentration:**
$$[H_0] = [H] + [HG] + 2[H_2] + 2[H_2G]$$

**Total guest concentration:**
$$[G_0] = [G] + [HG] + [H_2G]$$

### Algebraic Solution

From the equilibrium expressions:

$$[HG] = K_{HG}[H][G]$$
$$[H_2] = K_d[H]^2$$
$$[H_2G] = K_{H_2G}[H][HG] = K_{H_2G}K_{HG}[H]^2[G]$$

From guest mass balance:

$$[G_0] = [G] + K_{HG}[H][G] + K_{H_2G}K_{HG}[H]^2[G]$$
$$[G_0] = [G](1 + K_{HG}[H] + K_{H_2G}K_{HG}[H]^2)$$

Therefore:

$$[G] = \frac{[G_0]}{1 + K_{HG}[H] + K_{H_2G}K_{HG}[H]^2}$$

Substituting back:

$$[HG] = \frac{K_{HG}[H][G_0]}{1 + K_{HG}[H] + K_{H_2G}K_{HG}[H]^2}$$
$$[H_2G] = \frac{K_{H_2G}K_{HG}[H]^2[G_0]}{1 + K_{HG}[H] + K_{H_2G}K_{HG}[H]^2}$$

From host mass balance:

$$[H_0] = [H] + \frac{K_{HG}[H][G_0]}{1 + K_{HG}[H] + K_{H_2G}K_{HG}[H]^2} + 2K_d[H]^2 + 2\frac{K_{H_2G}K_{HG}[H]^2[G_0]}{1 + K_{HG}[H] + K_{H_2G}K_{HG}[H]^2}$$

Multiplying through by the denominator:

$$[H_0](1 + K_{HG}[H] + K_{H_2G}K_{HG}[H]^2) = [H](1 + K_{HG}[H] + K_{H_2G}K_{HG}[H]^2) + K_{HG}[H][G_0] + 2K_d[H]^2(1 + K_{HG}[H] + K_{H_2G}K_{HG}[H]^2) + 2K_{H_2G}K_{HG}[H]^2[G_0]$$

Expanding and collecting terms:

$$[H_0] + [H_0]K_{HG}[H] + [H_0]K_{H_2G}K_{HG}[H]^2 = [H] + K_{HG}[H]^2 + K_{H_2G}K_{HG}[H]^3 + K_{HG}[H][G_0] + 2K_d[H]^2 + 2K_dK_{HG}[H]^3 + 2K_dK_{H_2G}K_{HG}[H]^4 + 2K_{H_2G}K_{HG}[H]^2[G_0]$$

Rearranging into polynomial form:

$$2K_dK_{H_2G}K_{HG}[H]^4 + (K_{H_2G}K_{HG} + 2K_dK_{HG})[H]^3 + (2K_d + K_{HG} + 2K_{H_2G}K_{HG}[G_0] - K_{H_2G}K_{HG}[H_0])[H]^2 + (K_{HG}[G_0] - K_{HG}[H_0] + 1)[H] - [H_0] = 0$$

This leads to a **quartic polynomial equation** in [H] that typically requires **iterative numerical solution**.

### Iterative Solution Method

Due to the complexity of the algebraic solution, this system is typically solved iteratively:

1. **Initial guess**: Start with reasonable estimates for [H] and [G]
2. **Calculate complexes**: Use equilibrium expressions to find [HG], [H₂], [H₂G]
3. **Update free concentrations**: Use mass balance equations
4. **Check convergence**: Repeat until changes are below tolerance
5. **Calculate chemical shift**: Use final concentrations

### Chemical Shift Calculation (Guest-Observed)

For guest-observed NMR, only species containing guest molecules contribute:

$$\delta_{obs} = \frac{\delta_{G}[G] + \delta_{HG}[HG] + \delta_{H_2G}[H_2G]}{[G_0]}$$

The chemical shift change relative to free guest:

$$\Delta\delta = \delta_{obs} - \delta_G = \frac{\delta_{HG}[HG] + \delta_{H_2G}[H_2G]}{[G_0]}$$

### Physical Interpretation

This model describes systems where:
- **Host-guest binding** competes with **host dimerization**
- **Host dimers can also bind guest** molecules
- The system exhibits complex behavior depending on relative binding constants
- Higher host concentrations can favor either HG or H₂G formation depending on the equilibrium constants
- This model is particularly relevant for systems where host aggregation is significant

### Thermodynamic Relationships

The overall equilibrium for H₂G formation can be viewed as:

$$2H + G \rightleftharpoons H_2G$$

This can occur through two pathways:
1. **Sequential pathway**: $H + G \rightleftharpoons HG$, then $H + HG \rightleftharpoons H_2G$
2. **Alternative pathway**: $2H \rightleftharpoons H_2$, then $H_2 + G \rightleftharpoons H_2G$ (not present in this model)

The overall constant through the sequential pathway: $K_{overall} = K_{HG} \times K_{H_2G}$

This relationship helps in understanding the relative importance of the different pathways to complex formation. The formation of H₂G is directly dependent on the availability of both free host and the HG complex.

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

---

## References

1. **Fielding, L.** *Determination of association constants (Ka) from solution NMR data.* Tetrahedron 2000, 56, 6151-6170. [DOI](https://doi.org/10.1016/S0040-4020(00)00492-0)

2. **Hirose, K.** *A practical guide for the determination of binding constants.* J. Incl. Phenom. Macrocycl. Chem. 2001, 39, 193-209. [DOI](https://doi.org/10.1023/A:1011117412693)

3. **Connors, K.A.** *Binding Constants: The Measurement of Molecular Complex Stability.* John Wiley & Sons: New York, 1987. [DOI](https://doi.org/10.1002/bbpc.19870911223)

4. **Thordarson, P.** *Determining association constants from titration experiments in supramolecular chemistry.* Chem. Soc. 2011, 40. [DOI](https://doi.org/10.1039/C0CS00062K)

5. **Ulatowski, F.; et. al.** *Recognizing the Limited Applicability of Job Plots in Studying Host–Guest Interactions in Supramolecular Chemistry.* J. Org. Chem. 2016, 81, 1746-1756. [DOI](https://doi.org/10.1021/acs.joc.5b02909)
