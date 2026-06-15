# Reaction-Plane-Dependent Event Mixing

**A self-contained statistical method for extracting signals from noisy, correlated data** — developed for jet tomography of the quark-gluon plasma at RHIC, but the technique generalizes to any setting where the background structure is not independently measurable.

---

## The Problem

You want to measure correlations between particles emerging from a high-energy heavy ion collision. The signal of interest is an order of magnitude smaller than the background (from collective flow), and the background has the **same angular structure** as the signal. You can't fit it out. You can't measure it independently. Standard approaches rely on external reference data or ad-hoc assumptions.

**This framework measures and removes the entire background using only the data itself** — with no external inputs, no calibration runs, and no free parameters.

## The Method

### Step 1: Event Mixing — Isolate the Uncorrelated Background

Pair trigger particles with hadrons from *different* events in the same centrality class. Since the events are independent, flow correlations average out, leaving the **flat component**: $Y_\text{mix} = 4\pi \langle a_t\rangle\langle a_h\rangle$. But triggers favor higher-multiplicity events, so a correction $\kappa_a$ is needed.

### Step 2: Reaction-Plane Binned Mixing — Isolate the Correlated Background

Pair triggers with hadrons from events having the *same reaction-plane orientation*. This preserves the flow correlation: $Y_\text{mix}^{\psi} - Y_\text{mix} = 2\pi \langle\sigma b_t\rangle\langle\sigma b_h\rangle \cos(2\Delta\phi)$. But only the noisy *measured* reaction plane is available — $\kappa_b$ corrects for resolution damping.

### Step 3: In/Out-of-Plane Decomposition — Extract Geometry Dependence

Constrain triggers to align with or against the reaction plane. The difference isolates cross terms: $Y_\text{mix}^{\psi+} - Y_\text{mix}^{\psi-} = 8\langle a_h\rangle\langle b_t\sigma\rangle + 8\langle a_t\rangle\langle b_h\sigma\rangle \cos(2\Delta\phi)$. Two more corrections ($\kappa_{th}$, $\kappa_{ht}$) close the system.

### The Four $\kappa$ Corrections

| Factor | Corrects For | Formula |
|--------|-------------|---------|
| $\kappa_a$ | Multiplicity bias in flat term | $\langle a_t a_h\rangle / \langle a_t\rangle\langle a_h\rangle$ |
| $\kappa_b$ | Resolution damping in cosine term | $\langle b_t b_h\rangle \langle\sigma\rangle^2 / \langle b_t\sigma\rangle\langle b_h\sigma\rangle$ |
| $\kappa_{th}$ | Cross-term mismatch (t→h) | $\langle a_t b_h\sigma\rangle / \langle a_t\rangle\langle b_h\sigma\rangle$ |
| $\kappa_{ht}$ | Cross-term mismatch (h→t) | $\langle a_h b_t\sigma\rangle / \langle a_h\rangle\langle b_t\sigma\rangle$ |

Uncertainties propagated via 1000-trial Monte Carlo pseudo-experiments.

---

## Full Derivation

### 1. Single-Particle Ansatz

The azimuthal distribution of any particle species is determined by the second-order reaction plane $\psi_2$:

$$
f_i(\phi_i, \psi_{2i}) = a_i + b_i \cos\bigl(2(\phi_i - \psi_{2i})\bigr)
$$

$i \in \{\text{trigger } t,\ \text{hadron } h\}$. Coefficients $a_i$, $b_i$ depend on centrality ($N_\text{part}$) and $p_T$. Elliptic flow $v_{2i} = b_i / 2a_i$. Higher harmonics ($n>2$) are either analytically eliminated or negligible.

### 2. Same-Event Pair Background

The trigger–hadron $\Delta\phi = \phi_t - \phi_h$ distribution in the same event is the product of single-particle distributions with the same-event constraint $\delta(\psi_{2h}-\psi_{2t})$:

$$
F_{th}(\Delta\phi) = 4\pi\, a_t a_h + 2\pi\, b_t b_h \cos(2\Delta\phi)
$$

This is the **exact background**. First term isotropic; second oscillatory. Averaging over $N_\text{part}$:

$$
F_{th}(\Delta\phi) = 4\pi \langle a_t a_h\rangle + 2\pi \langle b_t b_h\rangle \cos(2\Delta\phi)
$$

### 3. Event Mixing — Measuring the Flat Term

Mix triggers from event $A$ with hadrons from event $B$ (same centrality, independent $\psi_2$). The cosine term vanishes:

$$
E_{th} = 4\pi \langle a_t\rangle\langle a_h\rangle
$$

This equals the flat term of $F_{th}$ **only if** $\langle a_t a_h\rangle = \langle a_t\rangle\langle a_h\rangle$. But the trigger rate increases with $N_\text{part}$, biasing triggered events toward larger $N_\text{part}$ than the average mixing event. The hadron yield also increases with $N_\text{part}$. Result: systematically underestimated flat background.

**Correction $\kappa_a$:**

$$
\kappa_a = \frac{\langle a_t a_h\rangle}{\langle a_t\rangle\langle a_h\rangle}
$$

Compute by estimating $a_t(N_\text{part})$ and $a_h(N_\text{part})$: plot per-centrality $\langle a_i\rangle$ vs $\langle N_\text{part}\rangle$, interpolate with cubic spline, integrate over $h(N_\text{part})$. Typical $\kappa_a = 1.01$–$1.05$.

After flat subtraction: $F_{th} - \kappa_a E_{th} = 2\pi \langle b_t b_h\rangle \cos(2\Delta\phi)$ — only the cosine term remains.

### 4. Reaction Plane Resolution

The measured reaction plane $\psi_{2}^{*}$ differs from the true $\psi_2$. The conditional distribution is a Gaussian periodic in $\pi$:

$$
g(\psi_2-\psi_{2}^{*}) = \frac{1}{\sqrt{2\pi\omega^2}} \sum_k \exp\!\left(-\frac{(\psi_2-\psi_{2}^{*}+k\pi)^2}{2\omega^2}\right)
$$

Convolving damps the flow coefficient by the resolution $\sigma = e^{-2\omega^2}$:

$$
f_{i}^{*} = \int g(\psi_2-\psi_{2}^{*})\, f_i\, d\psi_2 = a_i + \sigma b_i \cos\bigl(2(\phi_i - \psi_{2}^{*})\bigr)
$$

$\sigma$ is measured from the correlation between independent sub-detectors: $\sigma = \sqrt{\langle\cos(2(\psi_N-\psi_S))\rangle}$.

### 5. Reaction-Plane Binned Mixing — Measuring the Cosine Term

Mix triggers with hadrons whose measured $\psi_{2}^{*}$ is within $\varepsilon$ of the trigger's. Restricting $\psi_{2h}^{*}$ to $[\psi_{2t}^{*}-\varepsilon,\ \psi_{2t}^{*}+\varepsilon]$ preserves the flow correlation:

$$
E_{th}^{\psi} = 4\pi \langle a_t\rangle\langle a_h\rangle + 2\pi \langle\sigma b_t\rangle\langle\sigma b_h\rangle \cos(2\Delta\phi)
$$

(With $\varepsilon = \pi/80$, $\sin(2\varepsilon)/2\varepsilon \approx 1$.) Subtracting the unmixed result:

$$
E_{th}^{\psi} - E_{th} = 2\pi \langle\sigma b_t\rangle\langle\sigma b_h\rangle \cos(2\Delta\phi)
$$

But $\langle b_t b_h\rangle \neq \langle\sigma b_t\rangle\langle\sigma b_h\rangle / \langle\sigma\rangle^2$ — $\sigma$ correlates with $b$.

**Correction $\kappa_b$:**

$$
\kappa_b = \frac{\langle b_t b_h\rangle \, \langle\sigma\rangle^2}{\langle\sigma b_t\rangle\langle\sigma b_h\rangle}
$$

### 6. In/Out-of-Plane Backgrounds

Restrict triggers to those aligned with (in-plane) or perpendicular to (out-of-plane) the reaction plane:

$$
\psi_{2\pm}^{*}(\phi_t) = \phi_t + \alpha \pm \frac{\pi}{4}, \qquad \alpha \in \{0\ (\text{in}),\ \pi/2\ (\text{out})\}
$$

The same-event backgrounds:

$$
F_{th}^{\pm} = 2\pi\langle a_t a_h\rangle \pm 4\langle a_h b_t\sigma\rangle + \bigl(\pi\langle b_t b_h\rangle \pm 4\langle a_t b_h\sigma\rangle\bigr) \cos(2\Delta\phi)
$$

### 7. In/Out-of-Plane Event Mixing

Repeat RP-binned mixing with these limits:

$$
E_{th}^{\psi\pm} = 2\pi \langle a_t\rangle\langle a_h\rangle \pm 4\langle a_h\rangle\langle b_t\sigma\rangle \pm 4\langle a_t\rangle\langle b_h\sigma\rangle \cos(2\Delta\phi) + \pi \langle\sigma b_t\rangle\langle\sigma b_h\rangle \cos(2\Delta\phi)
$$

The difference isolates the cross terms:

$$
E_{th}^{\psi+} - E_{th}^{\psi-} = 8\langle a_h\rangle\langle b_t\sigma\rangle + 8\langle a_t\rangle\langle b_h\sigma\rangle \cos(2\Delta\phi)
$$

Two more corrections:

$$
\kappa_{th} = \frac{\langle a_t b_h\sigma\rangle}{\langle a_t\rangle\langle b_h\sigma\rangle}, \qquad
\kappa_{ht} = \frac{\langle a_h b_t\sigma\rangle}{\langle a_h\rangle\langle b_t\sigma\rangle}
$$

### 8. Corrected Backgrounds

All four $\kappa$ corrections applied:

**All-plane:**

$$
F_{th}^{*} = 4\pi\, \kappa_a \langle a_t\rangle\langle a_h\rangle + \frac{2\pi\, \kappa_b \langle\sigma b_t\rangle\langle\sigma b_h\rangle}{\langle\sigma\rangle^2} \cos(2\Delta\phi)
$$

**In/out-of-plane:**

$$
F_{th}^{\pm*} = 2\pi\, \kappa_a \langle a_t\rangle\langle a_h\rangle \pm 4\kappa_{ht} \langle a_h\rangle\langle b_t\sigma\rangle + \left(\frac{\pi\, \kappa_b \langle\sigma b_t\rangle\langle\sigma b_h\rangle}{\langle\sigma\rangle^2} \pm 4\kappa_{th} \langle a_t\rangle\langle b_h\sigma\rangle\right) \cos(2\Delta\phi)
$$

### 9. Final Subtraction

The per-trigger yield, corrected for detector acceptance $A(\Delta\phi)$:

$$
Y(\Delta\phi) = \frac{Y_\text{same}(\Delta\phi) - \kappa_a Y_\text{mix}(\Delta\phi) - \kappa_b \bigl(Y_\text{mix}^{\psi} - Y_\text{mix}\bigr) \cos(2\Delta\phi)}{A(\Delta\phi)}
$$

For in/out-of-plane, the $\pm$ cross terms are additionally removed using $\kappa_{th}$ and $\kappa_{ht}$.

### 10. Triangular Flow Elimination

Expanding to include $v_3$: $f_i = a_i + b_{2i}\cos(2\Delta\phi_i) + b_{3i}\cos(3\Delta\phi_i)$. The pair distribution acquires a $\cos(3\Delta\phi)$ term. By choosing away-side integration limits $|\Delta\phi - \pi| < \pi/3$, the $v_3$ contribution vanishes identically — no measurement, no model, no assumption.

---

## Subtraction in Action

### All-Plane

Raw same-event data (a) → after flat background subtraction (b) → final cleaned correlation (c):

<p align="center">
  <img src="figures/allplane_a.png" width="30%" />
  <img src="figures/allplane_b.png" width="30%" />
  <img src="figures/allplane_c.png" width="30%" />
</p>

### In-Plane

Five-step subtraction for triggers aligned with the reaction plane:

<p align="center">
  <img src="figures/inplane_a.png" width="18%" />
  <img src="figures/inplane_b.png" width="18%" />
  <img src="figures/inplane_c.png" width="18%" />
  <img src="figures/inplane_d.png" width="18%" />
  <img src="figures/inplane_e.png" width="18%" />
</p>

### Out-of-Plane

Five-step subtraction for triggers perpendicular to the reaction plane:

<p align="center">
  <img src="figures/outplane_a.png" width="18%" />
  <img src="figures/outplane_b.png" width="18%" />
  <img src="figures/outplane_c.png" width="18%" />
  <img src="figures/outplane_d.png" width="18%" />
  <img src="figures/outplane_e.png" width="18%" />
</p>

---

## Code Structure

```
src/                          C++ analysis framework
  gamAna.h/cc                 Event loop, track/photon handling
  gamAnaTree/                 Per-event tree building and histogramming
  EmcCluster.h/cc             EMCal cluster reconstruction and shower shapes
  CentralTrackData.h/cc       Track matching, DCA, particle ID
  TriggerParticle.h/cc        Trigger selection (photon, π⁰, asymmetric cuts)

macros/                       ~70 ROOT/C++ analysis and plotting macros
  artist_gamAnaTree.C         Main driver: Δφ histograms, event mixing
  acceptCorrect.C             Detector acceptance from RP-dependent mixing
  calcKComponents.C           κ corrections
  subtractKWeightedBG.C       Full background subtraction
  integAwaySide.C             Integrated yields and v₂ extraction
  CalcRPReso.C                Reaction plane resolution
  getFlow.C                   Flow background isolation
  draw*.C                     Publication figures
```

C++17 with ROOT. Processed ~20B minimum-bias Au+Au events over the grid (PROOF). Events are compressed offline for flexible cross-event mixing.

---

## About

**Z. Rowan**, "Removal of Anisotropic Background from Neutral Pion and Tagged Direct Photon–Hadron Correlations of Au+Au 200 GeV Collisions," PhD thesis, CUNY Graduate Center, 2022.

Advisor: Prof. Stefan Bathe. Built with the PHENIX Collaboration at RHIC, Brookhaven National Laboratory.
