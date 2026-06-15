# Photon-Hadron Azimuthal Correlations in Au+Au at $\sqrt{s_{NN}} = 200$ GeV

Reaction-plane-dependent event mixing for removing the collective background from two-particle correlations in heavy ion collisions. Measures the away-side hadron yield associated with $\pi^0$ and direct photon triggers to study path-length-dependent partonic energy loss in the quark-gluon plasma.

## Physics Context

Ultrarelativistic heavy ion collisions at RHIC create a QGP at energy densities $\gtrsim 1$ GeV/fm$^3$. Partons from initial hard scatterings lose energy via induced gluon bremsstrahlung (BDMPS/GLV) as they traverse this medium. The energy loss $\Delta E$ depends on the transport coefficient $\hat{q}$ and the path length $L$ through the medium: $\Delta E \propto \hat{q} L^2$ (BDMPS) or $\Delta E \propto \hat{q} L$ (GLV). Measuring this path-length dependence requires separating the signal from the overwhelming collective background.

## Background Removal: Full Derivation

### Single Particle Distribution

The azimuthal distribution of particles decomposes into an isotropic part plus elliptic flow relative to the reaction plane $\psi_2$:

$$f_i = a_i + b_i \cos\big(2(\phi_i - \psi_{2i})\big)$$

with the normalization constraint:

$$\int_0^{2\pi} d\phi_i\ \frac{1}{\pi} \int_0^{\pi} d\psi_{2i}\ f_i = 2\pi a_i$$

and the relation $b = 2 v_2 a$ connecting to the standard elliptic flow harmonic.

### Pair Background Convolution

The trigger-hadron pair distribution is the product of single-particle distributions with delta-function constraints enforcing same reaction plane and $\Delta\phi$ pairing:

$$\int_0^{2\pi} d\phi_t\ \frac{1}{\pi} \int_0^{\pi} d\psi_{2t}\ f_t \int_0^{2\pi} d\phi_h\ \delta(\Delta\phi - \phi_t + \phi_h) \int_0^{\pi} d\psi_{2h}\ \delta(\psi_{2h} - \psi_{2t})\ f_h$$

After full integration over the participant number distribution $h(N_{\text{part}})$, this yields:

$$F_{th} = 4\pi \langle a_t a_h \rangle + 2\pi \langle b_t b_h \rangle \cos(2\Delta\phi)$$

where $\langle x \rangle = \int_0^\infty dN_{\text{part}}\ h(N_{\text{part}})\ x$ is the expectation value over the participant distribution for a given centrality bin. The first term is the flat background; the second is the anisotropic $v_2$-modulated background. Both must be removed.

### Event Mixing Without Reaction Plane

Mixing triggers with hadrons from different events (same centrality, independent reaction planes) integrates out the cosine term:

$$E_{th} = 4\pi \langle a_t \rangle \langle a_h \rangle$$

This captures only the flat background, but with a bias: $\langle a_t \rangle \langle a_h \rangle \neq \langle a_t a_h \rangle$ because trigger particles preferentially emerge from higher-multiplicity events within each centrality bin. The correction factor is:

$$\kappa_a = \frac{\langle a_t a_h \rangle}{\langle a_t \rangle \langle a_h \rangle} = \frac{\int_0^\infty dN_{\text{part}}\ h(N_{\text{part}})\ a_t(N_{\text{part}})\ a_h(N_{\text{part}})}{\int_0^\infty dN_{\text{part}}\ h(N_{\text{part}})\ a_t(N_{\text{part}}) \int_0^\infty dN_{\text{part}}\ h(N_{\text{part}})\ a_h(N_{\text{part}})}$$

After flat subtraction: $F_{th} - \kappa_a E_{th} = 2\pi \langle b_t b_h \rangle \cos(2\Delta\phi)$.

### Reaction Plane Resolution

The measured reaction plane $\psi_2^*$ is related to the true plane $\psi_2$ via a periodic Gaussian resolution kernel:

$$g = \frac{1}{\sqrt{2\pi\omega^2}} \sum_{k=-\infty}^{\infty} \exp\left(-\frac{(\psi_2 - \psi_2^* + k\pi)^2}{2\omega^2}\right)$$

Convolving the single-particle distribution with this kernel damps the flow coefficient by the resolution:

$$f_i^* = a_i + \sigma\, b_i \cos\big(2(\phi_i - \psi_{2i}^*)\big), \qquad \sigma = e^{-2\omega^2}$$

The resolution is extracted from the correlation between independent sub-events (north/south detector halves):

$$\sigma = \langle \cos(2(\Psi^{\text{meas}} - \Psi^{\text{real}})) \rangle = \sqrt{\langle \cos(2(\Psi_N - \Psi_S)) \rangle}$$

### Event Mixing With Reaction Plane

Mixing triggers with hadrons from events having the same reaction plane orientation introduces a sinc-like factor from the finite bin width $\epsilon$:

$$E_{th}^{\psi} = 4\pi \langle a_t \rangle \langle a_h \rangle + 2\pi \frac{\sin(2\epsilon)}{2\epsilon} \langle \sigma b_t \rangle \langle \sigma b_h \rangle \cos(2\Delta\phi)$$

### In/Out-of-Plane Decomposition

Restricting triggers to those aligned in-plane or out-of-plane with the measured reaction plane:

$$F_{th}^{\pm} = 2\pi \langle a_t a_h \rangle \pm 4 \langle a_h b_t \sigma \rangle + \big(\pi \langle b_t b_h \rangle \pm 4 \langle a_t b_h \sigma \rangle\big) \cos(2\Delta\phi)$$

The cross terms are isolated from the difference of in and out-of-plane mixed event yields:

$$E_{th}^{p} - E_{th}^{m} = 8 \langle a_h \rangle \langle b_t \sigma \rangle + 8 \langle a_t \rangle \langle b_h \sigma \rangle \cos(2\Delta\phi)$$

### Complete $\kappa$ Correction Set

Four independent correction factors account for correlations between particle properties and the $N_{\text{part}}$ distribution:

$$\kappa_a = \frac{\langle a_t a_h \rangle}{\langle a_t \rangle \langle a_h \rangle}, \quad \kappa_b = \frac{\langle b_t b_h \rangle \langle \sigma \rangle^2}{\langle b_t \sigma \rangle \langle b_h \sigma \rangle}, \quad \kappa_{th} = \frac{\langle a_t b_h \sigma \rangle}{\langle a_t \rangle \langle b_h \sigma \rangle}, \quad \kappa_{ht} = \frac{\langle a_h b_t \sigma \rangle}{\langle a_h \rangle \langle b_t \sigma \rangle}$$

### Corrected Backgrounds

All-plane:

$$F_{th}^* = 4\pi \kappa_a \langle a_t \rangle \langle a_h \rangle + 2\pi \kappa_b \frac{\langle b_t \sigma \rangle \langle b_h \sigma \rangle}{\langle \sigma \rangle^2} \cos(2\Delta\phi)$$

In/out-of-plane:

$$F_{th}^{\pm*} = 2\pi \kappa_a \langle a_t \rangle \langle a_h \rangle \pm 4 \kappa_{ht} \langle a_h \rangle \langle b_t \sigma \rangle + \left(\pi \kappa_b \frac{\langle b_t \sigma \rangle \langle b_h \sigma \rangle}{\langle \sigma \rangle^2} \pm 4 \kappa_{th} \langle a_t \rangle \langle b_h \sigma \rangle\right) \cos(2\Delta\phi)$$

### Final Subtraction

The corrected per-trigger yield:

$$Y(\Delta\phi) = \frac{Y_{\text{same}}(\Delta\phi) - \kappa_a Y_{\text{mix}}(\Delta\phi) - \kappa_b (Y_{\text{mix}}^{\psi} - Y_{\text{mix}}) \cos(2\Delta\phi)}{A(\Delta\phi)}$$

where $A(\Delta\phi) = Y_{\text{mix}}(\Delta\phi) / n_{h,\text{exp}}$ is the acceptance function.

## Away-Side $v_2$ Extraction

### Generalized Flow Decomposition

The full two-particle correlation expanded to arbitrary flow harmonics:

$$n_{th} = 2\pi a_t a_h \left(1 + 2 \sum_{n=2}^{\infty} v_{tn} v_{hn} \cos(n\Delta\phi)\right)$$

### Exact $v_2$ Formula

The ratio of geometry-dependent yields isolates the associated particle $v_2$:

$$\frac{Y^+(\Delta\phi) - Y^-(\Delta\phi)}{Y(\Delta\phi)} = \frac{16\sigma\, v_{2a} \cos(2\Delta\phi) (a_t^2 - b_t^2)}{\pi (a_t^2 - 4\sigma^2 b_t^2/\pi^2)(1 + 4 v_{2a} b_t \cos(2\Delta\phi) / a_t)}$$

For small trigger $v_2$, this simplifies to:

$$\frac{\pi}{16\sigma} \frac{Y^+(\Delta\phi) - Y^-(\Delta\phi)}{Y(\Delta\phi)} \approx v_{2a} \cos(2\Delta\phi)$$

### Integrated $v_2$

From integrated yields over the away-side region $[2\pi/3, \pi]$:

$$\frac{\pi}{16\sigma} \frac{\int_{2\pi/3}^{\pi} Y^+(\Delta\phi) - \int_{2\pi/3}^{\pi} Y^-(\Delta\phi)}{(0.75 \pm 0.1) \int_{2\pi/3}^{\pi} Y(\Delta\phi)} = v_{2a}$$

The $0.75 \pm 0.1$ factor accounts for the shape dependence of the yield ratio.

## Direct Photon Tagging

### Tagging Matrix

The measured tagged photon counts form a linear system mixing true direct and decay contributions:

$$\begin{pmatrix} N_{\text{dec}}^t \\ N_{\text{dir}}^t \end{pmatrix} = \begin{pmatrix} \tau & \mu \\ 1-\mu & 1-\tau \end{pmatrix} \begin{pmatrix} N_{\text{dec}} \\ N_{\text{dir}} \end{pmatrix}$$

where $\tau$ is the tagging efficiency and $\mu$ is the mistag fraction. Inverting:

$$N_{\text{dir}} = \frac{\tau}{\tau - \mu} N_{\text{dir}}^t - \frac{1-\tau}{\tau - \mu} N_{\text{dec}}^t$$

### Decay Photon Spectrum

The decay photon energy distribution from a power-law pion spectrum $dN_\pi/dp_\pi \propto p_\pi^{-n}$:

$$\frac{dN_\gamma}{dE_\gamma} = \int_{p_{\pi,\min}}^{\infty} dp_\pi \frac{dN_\pi}{dp_\pi} \frac{1}{p_\pi} \propto \left(E_\gamma - \frac{m_\pi^2}{4E_\gamma}\right)^{-n} (\cos\theta)^{-n+1}$$

The fraction lost by the asymmetry cut $\alpha_{\max}$:

$$a \approx 1 - \left(\frac{2}{1+\alpha_{\max}}\right)^{-n}$$

### Composite Tagging Efficiency

Weighted over $\pi^0$ (80%) and $\eta$ (20%) decay channels:

$$\tau(p_T) = \frac{0.8 \left(N_{\text{tag}}^{\pi^0}/N_{\text{true}}^{\pi^0}\right) + 0.2 \left(N_{\text{tag}}^{\eta}/N_{\text{true}}^{\eta}\right)}{0.8 \left(N_{\text{meas}}^{\pi^0}/N_{\text{true}}^{\pi^0}\right) + 0.2 \left(N_{\text{meas}}^{\eta}/N_{\text{true}}^{\eta}\right)}$$

### Direct-to-Inclusive Ratio

$$R_\gamma = 1 + \frac{N_{\text{dir}}}{N_{\text{dec}}}$$

### Hybrid Yield Extraction

The direct photon-hadron yield from inclusive minus decay:

$$Y_{\text{dir}} = \frac{1}{N_{\text{dir}}} \left(N_{\text{inc}} Y_{\text{inc}} - N_{\text{dec}} Y_{\text{dec}}\right)$$

## Event Plane Construction

### Q-Vector and Recentering

The event plane angle from the second harmonic Q-vector:

$$Q_x^{\text{obs}} = \sum_i \cos(2\phi_i), \quad Q_y^{\text{obs}} = \sum_i \sin(2\phi_i), \quad \Psi^{\text{obs}} = \frac{1}{2} \tan^{-1}\left(\frac{Q_y^{\text{obs}}}{Q_x^{\text{obs}}}\right)$$

Detector non-uniformity corrected via recentering (z-score normalization):

$$Q_x^{\text{rec}} = \frac{Q_x^{\text{obs}} - \langle Q_x^{\text{obs}} \rangle}{\sigma_x}, \quad Q_y^{\text{rec}} = \frac{Q_y^{\text{obs}} - \langle Q_y^{\text{obs}} \rangle}{\sigma_y}$$

Residual asymmetries removed via Fourier flattening:

$$A_k = -\frac{2}{k} \langle \sin(k\Psi^{\text{rec}}) \rangle, \quad B_k = \frac{2}{k} \langle \cos(k\Psi^{\text{rec}}) \rangle$$

## Results

### $\pi^0$-Hadron Correlations

Per-trigger $\pi^0$-h yields from 0--60% Au+Au at $\sqrt{s_{NN}} = 200$ GeV. Away-side yield increases with trigger $p_T$; shape broadens at higher $\xi = \ln(p_{Tt}/p_{Th})$.

<!-- TODO: Add figure — Per-trigger π⁰-h yields, 0-60%, all p_Tt and ξ bins -->

### Centrality Comparison

Central (0--20%) vs mid-central (20--60%): away-side yield suppressed in central events, consistent with longer path length through the medium. Ratio approximately constant in $\xi$, with positive slope in trigger $p_T$.

<!-- TODO: Add figure — Central vs mid-central away-side comparison -->

### In vs Out-of-Plane

Positive away-side $v_2$ for all $\xi$ bins: in-plane hadrons (shorter path) enhanced relative to out-of-plane (longer path). Effect decreases at higher trigger $p_T$.

<!-- TODO: Add figure — In vs out-of-plane away-side comparison -->

### Direct Photon-Hadron

Central-to-mid-central ratio approximately constant as a function of trigger $p_T$ and $\xi$. Away-side $v_2$ hints at nonzero value.

<!-- TODO: Add figure — Direct photon-hadron results -->

## Code Structure

```
src/                          C++ analysis framework
  gamAna.h/cc                 SubsysReco module: event loop, track/photon handling
  gamAnaTree/                 Per-event tree building and histogram filling
  EmcCluster.h/cc             EMCal cluster reconstruction and shower shape cuts
  CentralTrackData.h/cc       CNT track matching, DCA, and PID (e/p, mass²)
  TriggerParticle.h/cc        Trigger selection (photon, π⁰, asymmetric cut)

macros/                       ROOT/C++ analysis and plotting
  artist_gamAnaTree.C         Main driver: builds Δφ histograms, applies mixing
  acceptCorrect.C             Acceptance function from RP-dependent mixing
  calcKComponents.C           κ_a, κ_b, κ_th, κ_ht computation
  subtractKWeightedBG.C       Full background subtraction
  integAwaySide.C             Integrated yield extraction and v₂
  CalcRPReso.C                Reaction plane resolution from sub-events
  getFlow.C                   Flow background isolation
  draw*.C                     Publication-quality figures
```

## References

- PHENIX Collaboration, PPG210 (direct photons in Au+Au)
- Sickles et al., Phys. Lett. B 693 (2010) 112
- PHENIX Collaboration, PPG023 (charged multiplicity)
