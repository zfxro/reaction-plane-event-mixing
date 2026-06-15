# Photon-Hadron Azimuthal Correlations in Au+Au Collisions at √s_NN = 200 GeV

Measurement of the away-side hadron yield associated with neutral pion and direct photon triggers in Au+Au collisions at RHIC/PHENIX, with reaction-plane-dependent background removal.

## Problem

In ultrarelativistic heavy ion collisions, a quark-gluon plasma (QGP) is created at energy densities exceeding 1 GeV/fm³. Partons produced in initial hard scatterings interact strongly with this QCD-colored medium, losing energy via gluon bremsstrahlung. This jet quenching is studied by measuring high-p_T hadrons correlated to trigger particles.

The largest obstacle is the overwhelming and partially correlated background from the collective behavior of the medium. Typically this background is subtracted with the help of results from an external analysis, or by assuming a zero yield at minimum (ZYAM). Both approaches have limitations:
- External flow measurements introduce systematic uncertainties from a separate analysis
- ZYAM provides no information about the reaction plane dependence of the signal

## Method: Reaction-Plane-Dependent Event Mixing

A technique is developed to remove the collective background from two-particle correlation measurements without relying on external flow measurements. The method exploits the reaction plane degree of freedom available in heavy ion collisions.

### Single Particle Distribution

The azimuthal distribution of particles is written as:

```
f_i = a_i + b_i cos(2(φ_i - ψ_2i))
```

where a_i and b_i are functions of centrality and p_T, φ_i is the particle azimuth, and ψ_2i is the second-order reaction plane angle of the event from which the particle emerged.

### Pair Background

The trigger-hadron pair azimuthal difference distribution, derived from the product of single particle distributions integrated over all angles, is:

```
F_th = 4π⟨a_t a_h⟩ + 2π⟨b_t b_h⟩ cos(2Δφ)
```

The first term is the flat (uncorrelated) background. The second term is the anisotropic (v_2-modulated) background. Both must be removed.

### Event Mixing Without Reaction Plane

Mixing triggers with hadrons from different events (same centrality) measures:

```
E_th = 4π⟨a_t⟩⟨a_h⟩
```

This captures only the flat background, scaled by a correction factor κ_a = ⟨a_t a_h⟩ / (⟨a_t⟩⟨a_h⟩) that accounts for the skewing of trigger events toward higher N_part within each centrality bin.

### Event Mixing With Reaction Plane

Mixing triggers with hadrons from events having the same reaction plane orientation measures:

```
E_th^ψ = 4π⟨a_t⟩⟨a_h⟩ + 2π⟨σ b_t⟩⟨σ b_h⟩ cos(2Δφ)
```

where σ = exp(-2ω²) is the reaction plane resolution. The difference isolates the cosine term:

```
E_th^ψ - E_th = 2π⟨σ b_t⟩⟨σ b_h⟩ cos(2Δφ)
```

### In and Out-of-Plane Backgrounds

By restricting triggers to those aligned in-plane or out-of-plane with respect to the measured reaction plane, the backgrounds become:

```
F_th^± = 2π⟨a_t a_h⟩ ± 4⟨a_h b_t σ⟩ + (π⟨b_t b_h⟩ ± 4⟨a_t b_h σ⟩) cos(2Δφ)
```

The cross terms are isolated from the difference of in and out-of-plane mixed event yields, with correction factors κ_th and κ_ht applied to account for expectation value discrepancies.

### Full Subtraction

The corrected backgrounds are subtracted from the same-event per-trigger yields:

```
Y(Δφ) = [Y_same(Δφ) - κ_a · Y_mix(Δφ) - κ_b · (Y_mix^ψ - Y_mix) · cos(2Δφ)] / A(Δφ)
```

where A(Δφ) = Y_mix(Δφ) / n_{h,exp} is the acceptance function computed from the mixed event yield divided by the expected hadron level.

## Results

### π⁰-Hadron Correlations

Per-trigger π⁰-h yields from 0-60% Au+Au collisions at √s_NN = 200 GeV show clear near-side (Δφ ≈ 0) and away-side (Δφ ≈ π) peaks. The away-side yield increases with trigger p_T while its shape broadens at higher ξ = ln(p_Tt/p_Th), consistent with increased fragmentation at lower associated momenta.

<!-- TODO: Add figure — Per-trigger π⁰-h yields, 0-60%, all p_Tt and ξ bins -->
<!-- Source: dphiplots/bigcan0_060_pi0_adj0.pdf, bigcan1_060_pi0_adj0.pdf -->

### Centrality Comparison: Path Length Dependence

Comparing central (0-20%) to mid-central (20-60%) collisions reveals a significant suppression of the away-side hadron yield in central events. Hadrons in central collisions traverse a larger average path through the medium, losing more energy. The ratio of central to mid-central integrated away-side yields is approximately constant as a function of ξ, with a positive slope as a function of trigger p_T suggesting that higher-momentum partons are less affected by the medium.

<!-- TODO: Add figure — Central vs mid-central away-side comparison -->
<!-- Source: dphiplots/bigcan_zoom_central_pi0_adj0.pdf -->

### In vs Out-of-Plane: Reaction Plane Dependence

A positive away-side v_2 is measured for all ξ bins, indicating that in-plane hadrons (shorter path through the medium) are enhanced relative to out-of-plane hadrons (longer path). The effect decreases at higher trigger p_T, consistent with the centrality comparison. Together, these results overwhelmingly suggest a path-length-dependent partonic energy loss in the medium.

<!-- TODO: Add figure — In vs out-of-plane away-side comparison -->
<!-- Source: dphiplots/bigcan_zoom_inout_pi0_adj0_cent2.pdf -->

### Direct Photon Tagging

A tagging method is introduced to measure the direct photon yield by removing decay photon contributions from the inclusive photon sample. The ratio R_γ of direct to inclusive photon yields is measured for various collision criteria and is consistent with previous PHENIX measurements. A hybrid method combining the plentiful inclusive photon sample with the tagging ratio produces statistically significant direct photon-hadron yields.

### Direct Photon-Hadron Results

The ratio of central to mid-central integrated away-side yields for direct photon triggers reveals an approximately constant suppression as a function of trigger p_T and associated hadron ξ. This is in contrast to previous I_AA measurements that showed an enhancement at lower momenta, though nuclear effects may cancel in the central-to-mid-central ratio. The away-side v_2 measurement hints at a nonzero value but with relatively large errors.

<!-- TODO: Add figure — Direct photon-hadron integrated away-side yield comparison -->
<!-- Source: integAwaySide or drawIntegAwaySide output -->

## Code Structure

```
src/
├── gamAna.h/cc              Main analysis framework (SubsysReco)
├── gamAnaTree/              Event-level tree building
│   ├── gamAnaTree.h         Analysis binning and histogram structure
│   └── gamAnaTree.cc        Event processing pipeline
├── EmcCluster.h/cc          EMCal cluster reconstruction
├── CentralTrackData.h/cc    CNT track matching and PID
├── ChargedTrack.h/cc        Base charged particle track
├── CentralArmTrack.h/cc     Central arm track with detector matching
├── TriggerParticle.h/cc     Trigger particle (photon/pi0 selection)
├── TrackCollection.h/cc     Templated track container
├── TrackCollectionBase.h/cc Abstract track container
└── gamConst.h/cc            Analysis constants and detector geometry

macros/
├── artist_gamAnaTree.C      Correlation analysis driver
├── acceptCorrect.C          Acceptance correction via event mixing
├── calcKComponents.C        K-factor background decomposition
├── subtractKWeightedBG.C    Weighted background subtraction
├── integAwaySide.C          Integrated away-side yield extraction
├── CalcRPReso.C             Reaction plane resolution
├── getFlow.C                Flow background extraction
├── drawDphiPlots.C          Δφ correlation plots
├── drawIntegAwaySide.C      Away-side yield figures
├── drawKComponents.C        K-component visualization
└── ...                      Additional analysis and plotting utilities
```

## References

- PHENIX Collaboration, "Direct photon production in Au+Au collisions at √s_NN = 200 GeV" (PPG210)
- A. Sickles et al., "Enhanced production of low p_T particles in di-jet events from Au+Au collisions at RHIC", Phys. Lett. B 693 (2010) 112-117
- PHENIX Collaboration, "Charged hadron multiplicity distributions in Au+Au collisions at RHIC" (PPG023)
