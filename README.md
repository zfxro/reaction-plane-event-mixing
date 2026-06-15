# Reaction-Plane-Dependent Event Mixing

**A self-contained statistical method for extracting signals from noisy, correlated data** — developed for jet tomography of the quark-gluon plasma at RHIC, but the technique generalizes to any setting where the background structure is not independently measurable.

[The Problem](#the-problem) | [The Method](#the-method) | [Full Derivation](#full-derivation) | [Code](#code-structure) | [Thesis](#thesis)

---

## The Problem

You want to measure correlations between particles emerging from a high-energy heavy ion collision. The signal you care about is an order of magnitude smaller than the background (from collective flow), and the background has the **same angular structure** as the signal. You can't fit it out. You can't measure it independently. Standard approaches rely on external reference data or ad-hoc assumptions.

**This framework measures and removes the entire background using only the data itself** — with no external inputs, no calibration runs, and no free parameters.

## The Method

### Step 1: Event Mixing — Isolate the Uncorrelated Background

Pair trigger particles with hadrons from *different* events in the same centrality class. Since the events are independent, flow correlations average out, leaving the **flat component**: `Y_mix = 4π ⟨a_t⟩⟨a_h⟩`. But triggers favor higher-multiplicity events, so a correction κₐ is needed.

### Step 2: Reaction-Plane Binned Mixing — Isolate the Correlated Background

Pair triggers with hadrons from events having the *same reaction-plane orientation*. This preserves the flow correlation: `Y_mix^ψ − Y_mix = 2π ⟨σ b_t⟩⟨σ b_h⟩ cos(2Δφ)`. But only the noisy *measured* reaction plane is available — κ_b corrects for resolution damping.

### Step 3: In/Out-of-Plane Decomposition — Extract Geometry Dependence

Constrain triggers to align with or against the reaction plane. The difference isolates cross terms: `Y_mix^ψ+ − Y_mix^ψ− = 8 ⟨a_h⟩⟨b_t σ⟩ + 8 ⟨a_t⟩⟨b_h σ⟩ cos(2Δφ)`. Two more corrections (κ_th, κ_ht) close the system.

<p align="center">
  <img src="figures/dphi_pi0_ptt0_sumcent1_dpsi0_xi0_a.png" width="32%" />
  <img src="figures/dphi_pi0_ptt0_sumcent1_dpsi0_xi0_b.png" width="32%" />
  <img src="figures/dphi_pi0_ptt0_sumcent1_dpsi0_xi0_c.png" width="32%" />
</p>
**Figure:** Background subtraction in action (π⁰–h, 20–60% centrality). Left: raw same-event data with flat background overlaid. Center: after flat subtraction, cosine background overlaid. Right: final cleaned correlation.

### The Four κ Corrections

| Factor | Corrects For | Formula | Typical Size |
|--------|-------------|---------|-------------|
| κₐ | Multiplicity bias in flat term | ⟨a_t a_h⟩ / ⟨a_t⟩⟨a_h⟩ | 1–5% |
| κ_b | Resolution damping in cosine term | ⟨b_t b_h⟩⟨σ⟩² / ⟨b_t σ⟩⟨b_h σ⟩ | 1–3% |
| κ_th | Cross-term mismatch (t→h) | ⟨a_t b_h σ⟩ / ⟨a_t⟩⟨b_h σ⟩ | ~2% |
| κ_ht | Cross-term mismatch (h→t) | ⟨a_h b_t σ⟩ / ⟨a_h⟩⟨b_t σ⟩ | ~2% |

Uncertainties propagated via 1000-trial Monte Carlo pseudo-experiments.

---

## Full Derivation

### 1. Single-Particle Ansatz

The azimuthal distribution of any particle species is determined by the second-order reaction plane ψ₂:

```
f_i(φ_i, ψ_{2i}) = a_i + b_i cos(2(φ_i − ψ_{2i}))
```

i ∈ {trigger t, hadron h}. Coefficients a_i, b_i depend on centrality (N_part) and p_T. Elliptic flow v_{2i} = b_i / 2a_i. Higher harmonics (n > 2) are either analytically eliminated or negligible.

Integrating over all angles gives the mean count per event:

```
∫ dφ_i · (1/π) ∫ dψ_{2i} f_i = n_i = 2π a_i
```

### 2. Same-Event Pair Background

The trigger–hadron Δφ = φ_t − φ_h distribution in the same event is the product of single-particle distributions with the same-event constraint δ(ψ_{2h}−ψ_{2t}):

```
F_th(Δφ) = ∫ dφ_t ∫ (dψ_{2t}/π) f_t · ∫ dφ_h ∫ dψ_{2h} δ(ψ_{2h}−ψ_{2t}) δ(Δφ−φ_t+φ_h) f_h
         = 4π a_t a_h + 2π b_t b_h cos(2Δφ)
```

This is the **exact background**. First term isotropic; second oscillatory. Averaging over N_part:

```
F_th(Δφ) = 4π ⟨a_t a_h⟩ + 2π ⟨b_t b_h⟩ cos(2Δφ)
```

### 3. Event Mixing — Measuring the Flat Term

Mix triggers from event A with hadrons from event B (same centrality, independent ψ₂). The cosine term vanishes:

```
E_th = 4π ⟨a_t⟩⟨a_h⟩
```

This equals the flat term of F_th **only if** ⟨a_t a_h⟩ = ⟨a_t⟩⟨a_h⟩. But the trigger rate increases with N_part, biasing triggered events toward larger N_part than the average mixing event. The hadron yield also increases with N_part. Result: systematically underestimated flat background.

**Correction κₐ:**

```
κ_a = ⟨a_t a_h⟩ / ⟨a_t⟩⟨a_h⟩
```

Compute by estimating a_t(N_part) and a_h(N_part): plot per-centrality ⟨a_i⟩ vs ⟨N_part⟩, interpolate with cubic spline, integrate over h(N_part). Typical κₐ = 1.01–1.05, largest in peripheral collisions where N_part distributions are most asymmetric.

After flat subtraction: `F_th − κ_a E_th = 2π ⟨b_t b_h⟩ cos(2Δφ)` — only the cosine term remains.

### 4. Reaction Plane Resolution

The measured reaction plane ψ₂* differs from the true ψ₂. The conditional distribution is a Gaussian periodic in π:

```
g(ψ₂−ψ₂*) = (1/√{2πω²}) Σ_k exp(−(ψ₂−ψ₂*+kπ)² / 2ω²)
```

Convolving damps the flow coefficient by the resolution σ = exp(−2ω²):

```
f_i* = ∫ g(ψ₂−ψ₂*) f_i dψ₂ = a_i + σ b_i cos(2(φ_i − ψ₂*))
```

The resolution σ is measured from the correlation between independent sub-detectors (north/south FVTX): σ = √⟨cos(2(ψ_N−ψ_S))⟩.

### 5. Reaction-Plane Binned Mixing — Measuring the Cosine Term

Mix triggers with hadrons whose measured ψ₂* is within ε of the trigger's. Restricting ψ_{2h}* to [ψ_{2t}*−ε, ψ_{2t}*+ε] preserves the flow correlation:

```
E_th^ψ = 4π ⟨a_t⟩⟨a_h⟩ + 2π · (sin(2ε)/2ε) · ⟨σ b_t⟩⟨σ b_h⟩ cos(2Δφ)
```

With ε = π/80, sin(2ε)/2ε ≈ 1. Subtracting the unmixed result:

```
E_th^ψ − E_th = 2π ⟨σ b_t⟩⟨σ b_h⟩ cos(2Δφ)
```

But ⟨b_t b_h⟩ ≠ ⟨σ b_t⟩⟨σ b_h⟩ / ⟨σ⟩² — σ correlates with b, so the product of averages differs from the average of the product.

**Correction κ_b:**

```
κ_b = ⟨b_t b_h⟩ ⟨σ⟩² / ⟨σ b_t⟩⟨σ b_h⟩
```

⟨σ b_t⟩ and ⟨σ b_h⟩ are extracted from E_th^ψ and the in/out-of-plane cross terms. b_i(N_part) is estimated by ⟨b_i⟩ ≈ ⟨σ b_i⟩/⟨σ⟩ (valid because σ varies weakly within a centrality bin), then spline-interpolated as in κₐ.

### 6. In/Out-of-Plane Backgrounds

Restrict triggers to those aligned with (in-plane) or perpendicular to (out-of-plane) the reaction plane:

```
ψ_{2±}*(φ_t) = φ_t + α ± π/4,   α ∈ {0 (in), π/2 (out)}
```

The same-event backgrounds with these limits:

```
F_th^± = 2π ⟨a_t a_h⟩ ± 4 ⟨a_h b_t σ⟩
       + (π ⟨b_t b_h⟩ ± 4 ⟨a_t b_h σ⟩) cos(2Δφ)
```

### 7. In/Out-of-Plane Event Mixing

Repeat RP-binned mixing with these limits:

```
E_th^{ψ±} = 2π ⟨a_t⟩⟨a_h⟩ ± 4 ⟨a_h⟩⟨b_t σ⟩ ± 4 ⟨a_t⟩⟨b_h σ⟩ cos(2Δφ)
          + π ⟨σ b_t⟩⟨σ b_h⟩ cos(2Δφ)
```

The difference isolates the cross terms:

```
E_th^{ψ+} − E_th^{ψ−} = 8 ⟨a_h⟩⟨b_t σ⟩ + 8 ⟨a_t⟩⟨b_h σ⟩ cos(2Δφ)
```

Two more corrections:

```
κ_{th} = ⟨a_t b_h σ⟩ / ⟨a_t⟩⟨b_h σ⟩,   κ_{ht} = ⟨a_h b_t σ⟩ / ⟨a_h⟩⟨b_t σ⟩
```

### 8. Measured Background Components

The four background terms are extracted directly from the mixing data:

<p align="center">
  <img src="figures/agah4pi_adj0_pi0_ptt0.png" width="45%" />
  <img src="figures/ahbgs8_adj0_pi0_ptt0.png" width="45%" />
</p>
<p align="center">
  <img src="figures/bgsbhs2pi_adj0_pi0_ptt0.png" width="45%" />
  <img src="figures/agbhs8_adj0_pi0_ptt0.png" width="45%" />
</p>
**Figure:** Measured background components for π⁰ (5 < p_T < 7 GeV) − h. Top left: flat term 4π⟨a_t a_h⟩. Top right: cross term 8⟨a_h b_t σ⟩. Bottom left: cosine term 2π⟨b_t b_h σ²⟩. Bottom right: cross term 8⟨a_t b_h σ⟩. Each panel shows the actual data measurement with fit.

The individual components ⟨a_t⟩, ⟨a_h⟩, ⟨b_t σ⟩, ⟨b_h σ⟩ are extracted from these four measured terms and used to compute the κ corrections.

### 9. Corrected Backgrounds

All four κ corrections applied:

**All-plane:**

```
F_th* = 4π κ_a ⟨a_t⟩⟨a_h⟩ + 2π κ_b ⟨σ b_t⟩⟨σ b_h⟩ / ⟨σ⟩² · cos(2Δφ)
```

**In/out-of-plane:**

```
F_th^{±*} = 2π κ_a ⟨a_t⟩⟨a_h⟩ ± 4 κ_{ht} ⟨a_h⟩⟨b_t σ⟩
          + (π κ_b ⟨σ b_t⟩⟨σ b_h⟩ / ⟨σ⟩² ± 4 κ_{th} ⟨a_t⟩⟨b_h σ⟩) cos(2Δφ)
```

### 10. Final Subtraction

The per-trigger yield, corrected for detector acceptance A(Δφ):

```
Y(Δφ) = [Y_same(Δφ) − κ_a Y_mix(Δφ) − κ_b (Y_mix^ψ − Y_mix) cos(2Δφ)] / A(Δφ)
```

For in/out-of-plane, the ± cross terms are additionally removed using κ_{th} and κ_{ht}.

### 11. Triangular Flow Elimination

Expanding to include v₃: `f_i = a_i + b_2i cos(2Δφ_i) + b_3i cos(3Δφ_i)`. The pair distribution acquires a cos(3Δφ) term. By choosing away-side integration limits |Δφ − π| < π/3, the v₃ contribution vanishes identically — no measurement, no model, no assumption. Just the geometry of the integral.

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

C++17 with ROOT. Processed ~20B minimum-bias Au+Au events over the grid (PROOF). Events are compressed offline for flexible cross-event mixing across the full data-taking period (~4 months, 2014).

---

## About

**Z. Rowan**, "Removal of Anisotropic Background from Neutral Pion and Tagged Direct Photon–Hadron Correlations of Au+Au 200 GeV Collisions," PhD thesis, CUNY Graduate Center, 2022.

Advisor: Prof. Stefan Bathe. Built with the PHENIX Collaboration at RHIC, Brookhaven National Laboratory.
