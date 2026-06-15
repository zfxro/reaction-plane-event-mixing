# Photon-Hadron Azimuthal Correlations in Au+Au at $\sqrt{s_{NN}} = 200$ GeV

Reaction-plane-dependent event mixing for removing the collective background from two-particle correlations in heavy ion collisions. Measures the away-side hadron yield associated with $\pi^0$ and direct photon triggers to study path-length-dependent partonic energy loss in the quark-gluon plasma.

## Physics Context

A fundamental constituent of matter, the quark, is confined within hadrons and not observed in isolation. At sufficiently high energy densities ($\gtrsim 1$ GeV/fm$^3$) the individuality of hadrons is blurred and quarks become deconfined. This newly discovered state of matter, known as quark-gluon plasma (QGP), existed in the early universe and is created in ultrarelativistic heavy ion collisions today. The quark gluon mixture rapidly thermalizes ($\lesssim 1$ fm/c after the heavy ion collision) before undergoing a cooling expansion and subsequent hadronization.

A proposed signature of the QGP, high transverse momentum ($p_T \gtrsim 2$ GeV) hadron suppression, was detected at RHIC in the early 2000's. Partons produced in initial hard scatterings interact strongly with the QCD colored medium losing energy via gluon bremsstrahlung. A quenching in central Au+Au collisions of high $p_T$ hadron spectra as compared to scaled spectra p+p collisions had been discovered.

Direct photon triggering has been referred to as the "golden channel" for jet tomography of the QGP. The colorless photon escapes the QCD medium unscathed, making its momentum a good proxy for the away-side jet's, and without introducing any surface bias.

The largest obstacle when making these types of measurements is the overwhelming and partially correlated background. Typically the background associated with the collective behavior of the medium is subtracted with the help of results from an external analysis. The remaining background is then subtracted by assuming a zero yield at minimum (ZYAM) or by implementing an absolute minimization method based on the event multiplicity. We expand on the latter by including additional information on the collision geometry. This allows for the removal of collective background without being reliant on previous measurements. An added bonus is the then straightforward measurement of in and out-of-plane yields.

## Background Removal: Full Derivation

### Single Particle Distribution

For nonzero impact parameter the overlapping region of the ions in the transverse plane resembles an ellipse. This causes an anisotropy of the azimuthal distribution of particles in the detector known as elliptic flow. The event by event multiplicity is modulated by a $\cos(2\Delta\phi)$ factor where $\Delta\phi = \phi - \psi_2$ (the difference in azimuth of a particle and the second order reaction plane).

The azimuthal distribution of particles decomposes into an isotropic part plus elliptic flow relative to the reaction plane $\psi_2$:

$$f_i = a_i + b_i \cos\big(2(\phi_i - \psi_{2i})\big)$$

with the normalization constraint:

$$\int_0^{2\pi} d\phi_i\ \frac{1}{\pi} \int_0^{\pi} d\psi_{2i}\ f_i = 2\pi a_i$$

and the relation $b = 2 v_2 a$ connecting to the standard elliptic flow harmonic.

### Pair Background Convolution

The delta function $\delta(\psi_{2h} - \psi_{2t})$ is required since the triggers and hadrons originate from the same event. The delta function $\delta(\Delta\phi - \phi_t + \phi_h)$ allows us to write this distribution in terms of the pair azimuthal difference. The trigger-hadron pair distribution is the product of single-particle distributions with delta-function constraints enforcing same reaction plane and $\Delta\phi$ pairing:

$$\int_0^{2\pi} d\phi_t\ \frac{1}{\pi} \int_0^{\pi} d\psi_{2t}\ f_t \int_0^{2\pi} d\phi_h\ \delta(\Delta\phi - \phi_t + \phi_h) \int_0^{\pi} d\psi_{2h}\ \delta(\psi_{2h} - \psi_{2t})\ f_h$$

After full integration over the participant number distribution $h(N_{\text{part}})$, this yields:

$$F_{th} = 4\pi \langle a_t a_h \rangle + 2\pi \langle b_t b_h \rangle \cos(2\Delta\phi)$$

where $\langle x \rangle = \int_0^\infty dN_{\text{part}}\ h(N_{\text{part}})\ x$ is the expectation value over the participant distribution for a given centrality bin. The first term is the flat background; the second is the anisotropic $v_2$-modulated background. Both must be removed.

### Event Mixing Without Reaction Plane

An attempt to measure the background can be made by mixing photons from the triggered events with hadrons of other events. Mixing triggers with hadrons from different events (same centrality, independent reaction planes) integrates out the cosine term:

$$E_{th} = 4\pi \langle a_t \rangle \langle a_h \rangle$$

This is precisely the first term of the pair background and we are halfway to removing the background. Unfortunately, we were naive in our approach. In the event mixing, a trigger is paired with hadrons from events having the same centrality as the triggered event. We have said that $a$ and $b$ are functions of centrality but to be more exact the physical quantity that they depend on is $N_{part}$, the number of nucleons participating in the collision. For each centrality there is not a single value of $N_{part}$, but a distribution. Because the trigger count increases with $N_{part}$, the triggered events are skewed toward values larger than the average $N_{part}$ of the mixed events. The end result is an underestimated level of background. The correction factor is:

$$\kappa_a = \frac{\langle a_t a_h \rangle}{\langle a_t \rangle \langle a_h \rangle} = \frac{\int_0^\infty dN_{\text{part}}\ h(N_{\text{part}})\ a_t(N_{\text{part}})\ a_h(N_{\text{part}})}{\int_0^\infty dN_{\text{part}}\ h(N_{\text{part}})\ a_t(N_{\text{part}}) \int_0^\infty dN_{\text{part}}\ h(N_{\text{part}})\ a_h(N_{\text{part}})}$$

After flat subtraction: $F_{th} - \kappa_a E_{th} = 2\pi \langle b_t b_h \rangle \cos(2\Delta\phi)$.

### Reaction Plane Resolution

The $\cos(2\Delta\phi)$ term of the background comes from correlations of the triggers and hadrons to the reaction plane. This term can potentially be measured if the events with which we mix triggers have the same reaction plane orientation as the triggered event. Unfortunately we falsely assumed that we had knowledge of the true reaction plane orientations. Of course this is not possible, but we do know the measured reaction plane. The relationship between the true, $\psi_2$, and measured, $\psi_2^*$, reaction plane angles is assumed to be normal with a period of $\pi$:

$$g = \frac{1}{\sqrt{2\pi\omega^2}} \sum_{k=-\infty}^{\infty} \exp\left(-\frac{(\psi_2 - \psi_2^* + k\pi)^2}{2\omega^2}\right)$$

Convolving the single-particle distribution with this kernel damps the flow coefficient by the resolution:

$$f_i^* = a_i + \sigma\, b_i \cos\big(2(\phi_i - \psi_{2i}^*)\big), \qquad \sigma = e^{-2\omega^2}$$

The resolution is extracted from the correlation between independent sub-events (north/south detector halves):

$$\sigma = \langle \cos(2(\Psi^{\text{meas}} - \Psi^{\text{real}})) \rangle = \sqrt{\langle \cos(2(\Psi_N - \Psi_S)) \rangle}$$

### Event Mixing With Reaction Plane

Since we have already removed the scalar term from the pair background we can isolate the cosine term by taking the difference of the mixing results with and without reaction plane binning. Mixing triggers with hadrons from events having the same reaction plane orientation introduces a sinc-like factor from the finite bin width $\epsilon$:

$$E_{th}^{\psi} = 4\pi \langle a_t \rangle \langle a_h \rangle + 2\pi \frac{\sin(2\epsilon)}{2\epsilon} \langle \sigma b_t \rangle \langle \sigma b_h \rangle \cos(2\Delta\phi)$$

The factor $\langle \sigma b_t \rangle \langle \sigma b_h \rangle$ can be found by fitting this difference in distributions by a $\cos(2\Delta\phi)$.

### In/Out-of-Plane Decomposition

Thus far, a $\psi_2$ dependence has only been introduced into the event mixing, in an attempt to measure the cosine term of the pair background. However the $\psi_2$ dependence of the signal could also be considered. We will look at two cases of the pair yield measurement: in and out-of-plane. Before doing so, the background in each case must be understood.

Restricting triggers to those aligned in-plane or out-of-plane with the measured reaction plane:

$$F_{th}^{\pm} = 2\pi \langle a_t a_h \rangle \pm 4 \langle a_h b_t \sigma \rangle + \big(\pi \langle b_t b_h \rangle \pm 4 \langle a_t b_h \sigma \rangle\big) \cos(2\Delta\phi)$$

The cross terms are what are new here. We isolate them by exploiting the symmetry of in and out-of-plane mixed event yields:

$$E_{th}^{p} - E_{th}^{m} = 8 \langle a_h \rangle \langle b_t \sigma \rangle + 8 \langle a_t \rangle \langle b_h \sigma \rangle \cos(2\Delta\phi)$$

The first, flat term can be measured by integrating over all $\Delta\phi$. The second, cosine term can be found by removing the previously calculated flat term and fitting the result with a simple cosine function.

### Complete $\kappa$ Correction Set

To account for the discrepancies between measured and actual background terms, four independent correction factors are calculated:

$$\kappa_a = \frac{\langle a_t a_h \rangle}{\langle a_t \rangle \langle a_h \rangle}, \quad \kappa_b = \frac{\langle b_t b_h \rangle \langle \sigma \rangle^2}{\langle b_t \sigma \rangle \langle b_h \sigma \rangle}, \quad \kappa_{th} = \frac{\langle a_t b_h \sigma \rangle}{\langle a_t \rangle \langle b_h \sigma \rangle}, \quad \kappa_{ht} = \frac{\langle a_h b_t \sigma \rangle}{\langle a_h \rangle \langle b_t \sigma \rangle}$$

In order to calculate $\kappa_{ij}$, $a_i$ and $b_j \sigma$ must be known as functions of $N_{part}$. Both $a_t(N_{part})$ and $a_h(N_{part})$ were previously determined. The products $b_t(N_{part}) \sigma(N_{part})$ and $b_h(N_{part}) \sigma(N_{part})$ can similarly be estimated from their expectation values, $\langle b_t \sigma \rangle$ and $\langle b_h \sigma \rangle$, which are extracted from the in/out-of-plane mixed event difference.

### Corrected Backgrounds

With measured background estimates $E_{th}$, $E_{th}^{\psi}$, and $E_{th}^{\pm}$; and correction factors $\kappa_a$, $\kappa_b$, $\kappa_{th}$, and $\kappa_{ht}$; we can begin chipping away at the inherent backgrounds. The measured all-plane background, including corrections, is:

$$F_{th}^* = 4\pi \kappa_a \langle a_t \rangle \langle a_h \rangle + 2\pi \kappa_b \frac{\langle b_t \sigma \rangle \langle b_h \sigma \rangle}{\langle \sigma \rangle^2} \cos(2\Delta\phi)$$

The in/out-of-plane backgrounds:

$$F_{th}^{\pm*} = 2\pi \kappa_a \langle a_t \rangle \langle a_h \rangle \pm 4 \kappa_{ht} \langle a_h \rangle \langle b_t \sigma \rangle + \left(\pi \kappa_b \frac{\langle b_t \sigma \rangle \langle b_h \sigma \rangle}{\langle \sigma \rangle^2} \pm 4 \kappa_{th} \langle a_t \rangle \langle b_h \sigma \rangle\right) \cos(2\Delta\phi)$$

### Final Subtraction

The corrected per-trigger yield:

$$Y(\Delta\phi) = \frac{Y_{\text{same}}(\Delta\phi) - \kappa_a Y_{\text{mix}}(\Delta\phi) - \kappa_b (Y_{\text{mix}}^{\psi} - Y_{\text{mix}}) \cos(2\Delta\phi)}{A(\Delta\phi)}$$

where $A(\Delta\phi) = Y_{\text{mix}}(\Delta\phi) / n_{h,\text{exp}}$ is the acceptance function.

## Away-Side $v_2$ Extraction

### Generalized Flow Decomposition

There are higher order effects caused by fluctuations in the initial parton density of the collision which introduce $\cos\big(n(\phi-\psi_n)\big)$ dependencies with $n > 2$. The full two-particle correlation expanded to arbitrary flow harmonics:

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

Neutral pions are responsible for producing 80% of the decay photons from 200 GeV Au+Au collisions. Nearly 99% of these pions decay via $\pi^0 \rightarrow 2\gamma$, a relatively simple channel to detect. It is with the help of this channel that we tag photons as product of decays. To reduce the combinatorial background, while maintaining an analytical understanding of the tagging efficiency, an asymmetry limit is imposed: $\alpha_{max} = 0.8$.

### Tagging Matrix

So in order to measure the direct photon signal just two functions must be known: the efficiency and mistag probability of the decay photon tagging method. Both will be obtained via simulations. The measured tagged photon counts form a linear system mixing true direct and decay contributions:

$$\begin{pmatrix} N_{\text{dec}}^t \\ N_{\text{dir}}^t \end{pmatrix} = \begin{pmatrix} \tau & \mu \\ 1-\mu & 1-\tau \end{pmatrix} \begin{pmatrix} N_{\text{dec}} \\ N_{\text{dir}} \end{pmatrix}$$

where $\tau$ is the tagging efficiency and $\mu$ is the mistag fraction. Inverting:

$$N_{\text{dir}} = \frac{\tau}{\tau - \mu} N_{\text{dir}}^t - \frac{1-\tau}{\tau - \mu} N_{\text{dec}}^t$$

### Decay Photon Spectrum

The opening angle between the two photons of a $\pi^0$ decay can extend to $180^{\circ}$. As we move to relativistic speeds, a smaller opening angle becomes increasingly more likely. The decay photon energy distribution from a power-law pion spectrum $dN_\pi/dp_\pi \propto p_\pi^{-n}$:

$$\frac{dN_\gamma}{dE_\gamma} = \int_{p_{\pi,\min}}^{\infty} dp_\pi \frac{dN_\pi}{dp_\pi} \frac{1}{p_\pi} \propto \left(E_\gamma - \frac{m_\pi^2}{4E_\gamma}\right)^{-n} (\cos\theta)^{-n+1}$$

PHENIX has measured $n$ in Au+Au 200 GeV to be $7.1 \pm 0.1$. The fraction lost by the asymmetry cut $\alpha_{\max}$:

$$a \approx 1 - \left(\frac{2}{1+\alpha_{\max}}\right)^{-n}$$

Thus around half, $a \approx 0.53$, of the decay photons are missed by the tagging method due to the max asymmetry cut. Although this is the largest factor that contributes to the tagging efficiency, it is not the only. Conversion probabilities, detector inefficiencies, and cluster merging effects must also be taken into account.

### Composite Tagging Efficiency

$\pi^0$'s are simulated flat in $p_T$, $\phi$, and $\eta$ and forced to decay to $2\gamma$. 80% of decay photons are produced via this channel. The majority of the remaining 20% are produced via $\eta \rightarrow 2\gamma$. These decays are simulated and run through the same procedure. Using these decays as a proxy for the entire remainder the tagging efficiency becomes:

$$\tau(p_T) = \frac{0.8 \left(N_{\text{tag}}^{\pi^0}/N_{\text{true}}^{\pi^0}\right) + 0.2 \left(N_{\text{tag}}^{\eta}/N_{\text{true}}^{\eta}\right)}{0.8 \left(N_{\text{meas}}^{\pi^0}/N_{\text{true}}^{\pi^0}\right) + 0.2 \left(N_{\text{meas}}^{\eta}/N_{\text{true}}^{\eta}\right)}$$

### Direct-to-Inclusive Ratio

A routinely referenced photon measurement is the ratio of inclusive to decay photons:

$$R_\gamma = 1 + \frac{N_{\text{dir}}}{N_{\text{dec}}}$$

### Hybrid Yield Extraction

The pure tagging method suffers from statistical limitations. A hybrid approach was developed in which per-trigger $\gamma_{\text{dec}} - h$ yields were produced from $\pi^0 - h$ yields with assistance from simulations. The direct photon-hadron yield is obtained from inclusive minus decay:

$$Y_{\text{dir}} = \frac{1}{N_{\text{dir}}} \left(N_{\text{inc}} Y_{\text{inc}} - N_{\text{dec}} Y_{\text{dec}}\right)$$

This subtraction method yields a higher statistical precision than the pure tagging method. The ingredients of the pure tagging method, the tagged decay and tagged direct $\gamma - h$ yields, have been usurped by the higher statistic $\pi^0 - h$ and $\gamma_{\text{inc}} - h$ yields.

## Event Plane Construction

Measurement of the second order reaction plane angle $\psi_2$, which facilitates the event mixing and background removal methods of this analysis, is done by the FVTX. The forward vertex detector tracks charged particles in the pseudorapidity region $|1.2| < \eta < |2.2|$. The north and south detectors each have four annular semiconductor layers.

### Q-Vector and Recentering

The reaction plane angle is measured by summing over all clusters in the FVTX while weighting by their position. The event plane angle from the second harmonic Q-vector:

$$Q_x^{\text{obs}} = \sum_i \cos(2\phi_i), \quad Q_y^{\text{obs}} = \sum_i \sin(2\phi_i), \quad \Psi^{\text{obs}} = \frac{1}{2} \tan^{-1}\left(\frac{Q_y^{\text{obs}}}{Q_x^{\text{obs}}}\right)$$

Detector non-uniformity corrected via recentering (z-score normalization):

$$Q_x^{\text{rec}} = \frac{Q_x^{\text{obs}} - \langle Q_x^{\text{obs}} \rangle}{\sigma_x}, \quad Q_y^{\text{rec}} = \frac{Q_y^{\text{obs}} - \langle Q_y^{\text{obs}} \rangle}{\sigma_y}$$

Residual asymmetries removed via Fourier flattening:

$$A_k = -\frac{2}{k} \langle \sin(k\Psi^{\text{rec}}) \rangle, \quad B_k = \frac{2}{k} \langle \cos(k\Psi^{\text{rec}}) \rangle$$

## Results

### $\pi^0$-Hadron Correlations

The high yield of the light meson and its relatively simple decay kinatics make the $\pi^0$ an ideal first candidate for testing the background removal machinery. The bias of the measured high $p_T$ neutral pions to have emerged from the surface of the nuclear interaction region, as well as the increased likelihood of the away-side jet to be parented by a gluon, make $\pi^0$ tagged hadrons an interesting probe of the QGP.

Per-trigger $\pi^0$-h yields from 0--60% Au+Au at $\sqrt{s_{NN}} = 200$ GeV. The near side amplitude is larger than the away side's. One possible explanation is the $\pi^0$'s surface bias. A high momentum pion is less likely to have originated near the center of the nuclear overlap region where it would be more susceptible to energy loss. The near side associated hadrons therefore have most likely traversed a shorter distance in the medium than those on the away side and consequently lost less energy. Away-side yield increases with trigger $p_T$; shape broadens at higher $\xi = \ln(p_{Tt}/p_{Th})$ as partons from the initial collision fragment to ever less momentum.

<!-- TODO: Add figure — Per-trigger π⁰-h yields, 0-60%, all p_Tt and ξ bins -->

### Centrality Comparison

A ratio of the per-trigger $\pi^0 - h$ yield in central to mid-central Au+Au collisions reveals a significant suppression of the away-side hadron yield. Integrated away-side yields, when plotted as a function of $\xi$, appear to lie along an exponential curve that is common to all trigger transverse momentum. Evidently the distribution of the away-side particle's fractional momentum ($p_{Th}/p_{Tt}$) is largely independent of the trigger's total momentum ($p_{Tt}$).

The ratio of central (0--20%) to mid-central (20--60%) integrated away-side yields allows for the observation of possible medium effects. In all $\xi$ bins a significant suppression of the away-side yield in central events is observed. A linear fit suggests a near constant effect in $\xi$. Plotted as a function of $p_{Tt}$, a seemingly constant suppression is again observed. A linear fit reveals a positive slope of some significance. The away-side yields of higher momentum trigger particles seem to be slightly less affected by the medium.

<!-- TODO: Add figure — Central vs mid-central away-side comparison -->

### In vs Out-of-Plane

A measure of the away-side $v_2$ provides additional evidence of a path length dependent partonic energy loss. Positive away-side $v_2$ for all $\xi$ bins: in-plane hadrons (shorter path through the medium) enhanced relative to out-of-plane (longer path). Effect decreases at higher trigger $p_T$.

<!-- TODO: Add figure — In vs out-of-plane away-side comparison -->

### Direct Photon-Hadron

Direct photon - hadron correlations are explored. To obtain the direct photon signal, a tagging method which removes contributions from decays is implemented. Statistical limitations implored us to develop the hybrid method. By leaning on the plentiful samples of inclusive photons and neutral pions, we are able to uncover a significant suppression of the direct photon triggered hadron yield. The reaction plane dependence is examined, where a nonzero value of the away-side hadron $v_2$ is discovered.

The ratio of central to mid-central integrated away-side yields for direct photon triggers reveals an approximately constant suppression as a function of trigger $p_T$ and associated hadron $\xi$. This is in contrast to the recent $I_{AA}$ measurement which shows an enhancement at lower momenta. It must though be noted that nuclear effects, to which an $I_{AA}$ measurement might be sensitive, may cancel in our ratio of per-trigger yields from central and mid-central A+A collisions.

Interestingly the suppression experienced in central Au+Au collisions by hadrons tagged by direct photons and neutral pions is roughly the same. Naively one would expect a greater suppression of per-trigger $\pi^0 - h$ away side yields due to the pion's surface bias and away-side jet's tendency to be gluonic. A direct photon has no such bias, and its away-side jet is more likely to have been spawned by a quark, which interact less strongly with the quark-gluon plasma.

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
