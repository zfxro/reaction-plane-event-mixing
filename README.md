# Photon–Hadron Correlation Analysis Framework

A reaction-plane-dependent event mixing method that measures and removes the collective background from two-particle correlations without any external flow measurements — elliptic flow ($v_2$) is determined from the same data via sub-event resolution and a four-parameter $\kappa$ correction framework. Triangular flow ($v_3$) contributions are eliminated analytically through the choice of integration limits. An in/out-of-plane decomposition extracts the away-side hadron $v_2$, isolating the geometry dependence of partonic energy loss. A hybrid tagging method (inclusive $\gamma$ minus $\pi^0$-tagged decay $\gamma$) overcomes the statistical limitations of pure decay-photon tagging, enabling statistically significant direct photon-hadron correlation measurements.

**Contents**

- [Physics Context](#physics-context)
- [Experiment and Data](#experiment-and-data)
- [Background Removal: Full Derivation](#background-removal-full-derivation)
- [Neutral Pion - Hadron Correlations](#neutral-pion---hadron-correlations)
- [Direct Photon Tagging](#direct-photon-tagging)
- [Direct Photon - Hadron Correlations](#direct-photon---hadron-correlations)
- [Discussion](#discussion)
- [Code Structure](#code-structure)

---

## Physics Context

A fundamental constituent of matter, the quark, is confined within hadrons and not observed in isolation. At sufficiently high energy densities ($\gtrsim 1$ GeV/fm$^3$) the individuality of hadrons is blurred and quarks become deconfined. This newly discovered state of matter, known as quark-gluon plasma (QGP), existed in the early universe and is created in ultrarelativistic heavy ion collisions today. Microseconds after the Big Bang, the universe was a primordial soup of quarks and gluons with $T > 10^{12}$ K (the temperature corresponding to the fundamental energy scale of QCD, $\Lambda_{QCD} \sim 200$ MeV). These conditions can be recreated at the Relativistic Heavy Ion Collider (RHIC) at Brookhaven National Laboratory where measured energy densities of "little bangs" have exceeded 1 GeV/fm$^3$. The quark gluon mixture rapidly thermalizes ($\lesssim 1$ fm/c after the heavy ion collision) before undergoing a cooling expansion and subsequent hadronization.

A proposed signature of the QGP, high transverse momentum ($p_T \gtrsim 2$ GeV) hadron suppression, was detected at RHIC in the early 2000's. Partons produced in initial hard scatterings interact strongly with the QCD colored medium losing energy via gluon bremsstrahlung. A quenching in central Au+Au collisions of high $p_T$ hadron spectra as compared to scaled spectra p+p collisions had been discovered.

In addition to the $p_T$ spectra of singular particles, dihadron correlation measurements make a useful tool for probing the QGP. Triggering on high $p_T$ particles increases the fraction of measured hadrons belonging to collimated jets that have originated from initial hard scatterings. Comparisons of per-trigger yields from A+A and p+p collisions do not require a scaling by the number of individual nucleon collisions and have provided evidence of medium-induced partonic energy loss.

Direct photon triggering has been referred to as the "golden channel" for jet tomography of the QGP. The colorless photon escapes the QCD medium unscathed, making its momentum a good proxy for the away-side jet's, and without introducing any surface bias.

The largest obstacle when making these types of measurements is the overwhelming and partially correlated background. Typically the background associated with the collective behavior of the medium is subtracted with the help of results from an external analysis. The remaining background is then subtracted by assuming a zero yield at minimum (ZYAM) or by implementing an absolute minimization method based on the event multiplicity. We expand on the latter by including additional information on the collision geometry. This allows for the removal of collective background without being reliant on previous measurements. An added bonus is the then straightforward measurement of in and out-of-plane yields.

This background removal method is first implemented in the measurement of $\pi^0 - h$ correlations. The high yield of the light meson and its relatively simple decay kinematics ($\pi^0 \rightarrow 2\gamma$, branching ratio $98.823 \pm 0.034\%$) make it an ideal first candidate. The bias of the measured high $p_T$ neutral pions to have emerged from the surface of the nuclear interaction region, as well as the increased likelihood of the away-side jet to be parented by a gluon, make $\pi^0$ tagged hadrons an interesting probe of the QGP.

Direct photon - hadron correlations are then explored. To obtain the direct photon signal, a tagging method that removes contributions from decays is implemented. A hybrid approach using inclusive photons and neutral pions uncovers a significant suppression of the direct photon triggered hadron yield. The reaction plane dependence is examined, where a nonzero value of the away-side hadron $v_2$ is discovered.

---

## Experiment and Data

### PHENIX at RHIC

The analysis is conducted on the Run 14 (2014) Au+Au $\sqrt{s_{NN}} = 200$ GeV dataset, for which PHENIX recorded over 20 billion minimum bias events. The PHENIX detector consists of several subdetectors. Global detectors at forward rapidity characterize the collisions: the Beam-Beam Counters (BBC, $|3.0| < \eta < |3.9|$) measure collision centrality via the total charge deposited, and the Forward Vertex Detector (FVTX, $|1.2| < \eta < |2.2|$) measures the second-order reaction plane. Each of the two central arm spectrometers at mid-rapidity ($|\eta| < 0.35$) covers approximately one steradian in azimuth. The Drift Chamber (DC) provides charged particle momentum reconstruction. The Electromagnetic Calorimeters (EMCal, PbSc and PbGl) are used for photon measurement and triggering. The Pad Chambers (PC) provide charged particle tracking and, together with the Ring Imaging Cherenkov detector (RICH), enable electron/hadron separation.

### Collision Centrality

A pair of Beam-Beam Counters surround the beam pipe 1.4 m from the nominal collision vertex, covering $|3.0| < \eta < |3.9|$ over the full azimuth. A greater amount of charge is deposited into the BBC from head-on collisions than from glancing ones. Percentiles of the total charge distribution describe the centrality of the collision. Monte Carlo Glauber simulations relate the measured centrality to physical quantities such as the number of participating nucleons, $N_{part}$, and the impact parameter. The $N_{part}$ distributions for the centrality bins used in this analysis are obtained from these simulations.

### Reaction Plane Measurement

For nonzero impact parameter the overlapping region of the ions in the transverse plane resembles an ellipse, causing an anisotropy of the azimuthal distribution of particles known as elliptic flow. Measurement of the second order reaction plane angle $\psi_2$, which facilitates the event mixing and background removal methods, is done by the FVTX. The north and south detectors each have four annular semiconductor layers. The reaction plane angle is measured by summing over all clusters in the FVTX while weighting by their position:

$$Q_x^{\text{obs}} = \sum_i \cos(2\phi_i),\qquad Q_y^{\text{obs}} = \sum_i \sin(2\phi_i),\qquad \Psi^{\text{obs}} = \frac{1}{2} \tan^{-1}\left(\frac{Q_y^{\text{obs}}}{Q_x^{\text{obs}}}\right)$$

Detector non-uniformities are corrected via recentering (z-score normalization), and residual asymmetries are removed via Fourier flattening:

$$Q_x^{\text{rec}} = \frac{Q_x^{\text{obs}} - \langle Q_x^{\text{obs}} \rangle}{\sigma_x},\qquad Q_y^{\text{rec}} = \frac{Q_y^{\text{obs}} - \langle Q_y^{\text{obs}} \rangle}{\sigma_y}$$

$$A_k = -\frac{2}{k} \langle \sin(k\Psi^{\text{rec}}) \rangle,\qquad B_k = \frac{2}{k} \langle \cos(k\Psi^{\text{rec}}) \rangle$$

The reaction plane resolution is extracted from the correlation between independent sub-events (north/south detector halves):

$$\sigma = \langle \cos(2(\Psi^{\text{meas}} - \Psi^{\text{real}})) \rangle = \sqrt{\langle \cos(2(\Psi_N - \Psi_S)) \rangle}$$

### Charged Hadron Measurement

Charged hadrons are measured in the Drift Chamber at mid-rapidity ($|\eta| < 0.35$). The following track cuts are applied: $0.5 \leq p_T \leq 10$ GeV, DC track quality requirements, 3$\sigma$ EMCal and PC3 matching, RICH veto ($n0 < 0$ for $p_T < 5$ GeV) to remove electrons, and a $dep < -2.0$ cut on the sigmalized energy-over-momentum ratio to further reject electrons.

### Photon Measurement

Photons are measured in the Electromagnetic Calorimeters (PbSc and PbGl). Trigger photons are selected with $5.0 \leq p_T \leq 15.0$ GeV, a minimum cluster energy $ecore \geq 0.5$ GeV, a shower shape cut $prob\_photon \geq 0.02$, a PC3 charge veto for electron rejection, and a fiducial cut ensuring that decay partner photons land within the detector acceptance. Partner photons used for $\pi^0$ tagging are identified by a pair invariant mass within $3\sigma$ of the $\pi^0$ mass ($|M_{pair} - 0.138| \leq 0.03$ GeV) and an asymmetry cut $\alpha < 0.8$.

### Event Mixing Procedure

Mixed $\Delta\phi$ distributions are created by pairing trigger particles with charged hadrons of different events having the same:

- **run number**: mitigates detector response fluctuations over the data-taking period
- **centrality**: governs the multiplicity and overall background level; 5% and 10% wide bins are used
- **z-vertex**: $\pm 15$ cm, divided into 6 bins of 5 cm
- **reaction plane**: if enabled, the measured $\psi_2$ of the mixing event must be within $\pi/80$ of the trigger event's $\psi_2$, ensuring $\sin(2\epsilon)/2\epsilon \approx 1$

Events are compressed and saved for mixing offline to allow triggers to be mixed with events from the entire run. The mixed event per-trigger yield is:

$$\Upsilon_m(\Delta\phi) = \frac{1}{N_t} \sum_{i=1}^{N_t} \frac{1}{N_{mix,i}}\sum_{j=1}^{N_{mix,i}} N_{h,j}(\Delta\phi)$$

An acceptance function $A(\Delta\phi) = \Upsilon_m(\Delta\phi) / n_{h,exp}$ is constructed using the expected charged hadron yield from previously published PHENIX measurements (PPG023), and the acceptance-corrected yields are denoted $Y_{mix}(\Delta\phi)$.

---

## Background Removal: Full Derivation

### Single Particle Azimuthal Distribution

The azimuthal distributions of both trigger and background particles can be written in the form:

$$f_i = a_i + b_i \cos\big(2(\phi_i - \psi_{2i})\big)$$

where $i = t$ for triggers and $i = h$ for hadrons. $a_i$ and $b_i$ are functions of centrality and $p_T$ (or $\xi = \ln(p_{Tt}/p_{Th})$). While higher order reaction plane dependencies ($n > 2$) exist, the away-side integration limits ($|\Delta\phi - \pi| < \pi/3$) are chosen to eliminate $v_3$ contributions, and ever higher order terms are negligible.

Integrating $f_i$ over all $\phi_i$ and $\psi_i$ relates $a_i$ to the average number of particles per event: $n_i = 2\pi a_i$. The elliptic flow harmonic is related by $b_i = 2 v_{2i} a_i$.

### Pair Background Convolution

The trigger-hadron pair $\Delta\phi$ distribution is calculated by taking the product of single particle azimuthal distributions while integrating over all angles:

$$\int_0^{2\pi} d\phi_t\ \frac{1}{\pi} \int_0^{\pi} d\psi_{2t}\ f_t \int_0^{2\pi} d\phi_h\ \delta(\Delta\phi - \phi_t + \phi_h) \int_0^{\pi} d\psi_{2h}\ \delta(\psi_{2h} - \psi_{2t})\ f_h$$

The $\delta(\psi_{2h} - \psi_{2t})$ constrains triggers and hadrons to the same event; $\delta(\Delta\phi - \phi_t + \phi_h)$ writes the distribution in terms of the pair azimuthal difference. After integration:

$$F_{th} = 4\pi a_t a_h + 2\pi b_t b_h \cos(2\Delta\phi)$$

Folding over $\Delta\phi = \pi$ (since both signal and background are $\pi$-periodic) and restricting $\Delta\phi \in [0,\pi)$ gives:

$$4\pi a_t a_h + 2\pi b_t b_h \cos(2\Delta\phi)$$

This is the background that must be removed. By including the $N_{part}$ distribution $h(N_{part})$ and integrating, the full background becomes:

$$F_{th} = 4\pi \langle a_t a_h \rangle + 2\pi \langle b_t b_h \rangle \cos(2\Delta\phi)$$

where $\langle x \rangle = \int_0^\infty dN_{part}\ h(N_{part})\ x$.

### Event Mixing Without Reaction Plane

An attempt to measure the background can be made by mixing triggers with hadrons of other events. Mixing triggers with hadrons from different events (same centrality, independent reaction planes) integrates out the $\cos(2\Delta\phi)$ term:

$$E_{th} = 4\pi \langle a_t \rangle \langle a_h \rangle$$

This is precisely the first term of the pair background — we are halfway to removing it. However, because the trigger count increases with $N_{part}$, triggered events are skewed toward larger $N_{part}$ than the average of the mixed events. The end result is an underestimated background. A correction factor $\kappa_a$ is required:

$$\kappa_a = \frac{\langle a_t a_h \rangle}{\langle a_t \rangle \langle a_h \rangle}$$

After flat subtraction: $F_{th} - \kappa_a E_{th} = 2\pi \langle b_t b_h \rangle \cos(2\Delta\phi)$.

To calculate $\kappa_a$, $a_t(N_{part})$ and $a_h(N_{part})$ are estimated by plotting $\langle a \rangle$ as a function of $\langle N_{part} \rangle$ in each centrality bin and interpolating with a third order spline. The estimation works well in central collisions where the $N_{part}$ distributions are symmetric about their means. In peripheral events (centrality $> 60\%$) the asymmetry of the distributions causes a breakdown in the approximation, and these collisions are excluded from the analysis.

### Reaction Plane Resolution

The $\cos(2\Delta\phi)$ term of the background comes from correlations of the triggers and hadrons to the reaction plane. This term can be measured if the events with which we mix triggers have the same reaction plane orientation as the triggered event. We do not know the true reaction plane, but we know the measured one. The relationship between the true ($\psi_2$) and measured ($\psi_2^*$) reaction plane angles is assumed to be normal with period $\pi$:

$$g = \frac{1}{\sqrt{2\pi\omega^2}} \sum_{k=-\infty}^{\infty} \exp\left(-\frac{(\psi_2 - \psi_2^* + k\pi)^2}{2\omega^2}\right)$$

Convolving the single-particle distribution with this kernel damps the flow coefficient by the resolution:

$$f_i^* = a_i + \sigma\, b_i \cos\big(2(\phi_i - \psi_{2i}^*)\big), \qquad \sigma = e^{-2\omega^2}$$

### Event Mixing With Reaction Plane

Repeating the event mixing while requiring the measured reaction plane of mixing events to be within $\epsilon$ of the triggered event's:

$$E_{th}^{\psi} = 4\pi \langle a_t \rangle \langle a_h \rangle + 2\pi \frac{\sin(2\epsilon)}{2\epsilon} \langle \sigma b_t \rangle \langle \sigma b_h \rangle \cos(2\Delta\phi)$$

With sufficiently fine RP bins ($\epsilon = \pi/80$) such that $\sin(2\epsilon)/2\epsilon \approx 1$:

$$E_{th}^{\psi} = 4\pi \langle a_t \rangle \langle a_h \rangle + 2\pi \langle \sigma b_t \rangle \langle \sigma b_h \rangle \cos(2\Delta\phi)$$

Since the scalar term is already removed, the cosine term is isolated by taking the difference:

$$E_{th}^{\psi} - E_{th} = 2\pi \langle \sigma b_t \rangle \langle \sigma b_h \rangle \cos(2\Delta\phi)$$

A correction $\kappa_b$ is required for the discrepancy between measured and actual cosine terms:

$$\kappa_b = \frac{\langle b_t b_h \rangle \langle \sigma \rangle^2}{\langle \sigma b_t \rangle \langle \sigma b_h \rangle}$$

### In/Out-of-Plane Decomposition

The $\psi_2$ dependence of the *signal* can also be considered. Restricting triggers to those aligned in-plane or out-of-plane with the measured reaction plane (with suitable integration limits $\psi_{2\pm}^*(\phi_t) = \phi_t + \alpha \pm \pi/4$, where $\alpha = 0$ for in-plane and $\alpha = \pi/2$ for out-of-plane):

$$F_{th}^{\pm} = 2\pi \langle a_t a_h \rangle \pm 4 \langle a_h b_t \sigma \rangle + \big(\pi \langle b_t b_h \rangle \pm 4 \langle a_t b_h \sigma \rangle\big) \cos(2\Delta\phi)$$

The cross terms $\langle a_h b_t \sigma \rangle$ and $\langle a_t b_h \sigma \rangle$ are isolated by exploiting the symmetry of in and out-of-plane mixed event yields:

$$E_{th}^{p} - E_{th}^{m} = 8 \langle a_h \rangle \langle b_t \sigma \rangle + 8 \langle a_t \rangle \langle b_h \sigma \rangle \cos(2\Delta\phi)$$

The flat term is obtained by integrating over all $\Delta\phi$; the cosine term is found by removing the flat term and fitting with $\cos(2\Delta\phi)$. Two additional corrections $\kappa_{th}$ and $\kappa_{ht}$ account for discrepancies:

$$\kappa_{ij} = \frac{\langle a_i b_j \sigma \rangle}{\langle a_i \rangle \langle b_j \sigma \rangle}$$

### Complete $\kappa$ Correction Set

Four independent correction factors are calculated:

$$\kappa_a = \frac{\langle a_t a_h \rangle}{\langle a_t \rangle \langle a_h \rangle},\qquad \kappa_b = \frac{\langle b_t b_h \rangle \langle \sigma \rangle^2}{\langle b_t \sigma \rangle \langle b_h \sigma \rangle},\qquad \kappa_{th} = \frac{\langle a_t b_h \sigma \rangle}{\langle a_t \rangle \langle b_h \sigma \rangle},\qquad \kappa_{ht} = \frac{\langle a_h b_t \sigma \rangle}{\langle a_h \rangle \langle b_t \sigma \rangle}$$

All $\kappa$ corrections are on the order of a few percent. Monte Carlo simulations (1000 pseudo-experiments sampling the measured data points within their statistical errors) are conducted to produce the error on each $\kappa$ value.

### Corrected Backgrounds

The measured all-plane background, including corrections:

$$F_{th}^* = 4\pi \kappa_a \langle a_t \rangle \langle a_h \rangle + 2\pi \kappa_b \frac{\langle b_t \sigma \rangle \langle b_h \sigma \rangle}{\langle \sigma \rangle^2} \cos(2\Delta\phi)$$

The in/out-of-plane backgrounds:

$$F_{th}^{\pm*} = 2\pi \kappa_a \langle a_t \rangle \langle a_h \rangle \pm 4 \kappa_{ht} \langle a_h \rangle \langle b_t \sigma \rangle + \left(\pi \kappa_b \frac{\langle b_t \sigma \rangle \langle b_h \sigma \rangle}{\langle \sigma \rangle^2} \pm 4 \kappa_{th} \langle a_t \rangle \langle b_h \sigma \rangle\right) \cos(2\Delta\phi)$$

The analysis is conducted in 5% centrality bins (0--5, 5--10, ..., 50--60%) which are later merged into central (0--20%), mid-central (20--60%), and inclusive (0--60%) bins by weighting the $\kappa$-corrected backgrounds by the number of triggers in each fine centrality bin.

### Final Subtraction

The corrected per-trigger yield:

$$Y(\Delta\phi) = \frac{Y_{\text{same}}(\Delta\phi) - \kappa_a Y_{\text{mix}}(\Delta\phi) - \kappa_b (Y_{\text{mix}}^{\psi} - Y_{\text{mix}}) \cos(2\Delta\phi)}{A(\Delta\phi)}$$

where $A(\Delta\phi) = Y_{\text{mix}}(\Delta\phi) / n_{h,\text{exp}}$ is the acceptance function derived from previously published charged hadron yields.

---

## Neutral Pion - Hadron Correlations

### $\pi^0$ Decay Kinematics

The opening angle between the two photons of a $\pi^0$ decay can extend to $180^\circ$. The energy distribution of a decay photon in the lab frame is:

$$\frac{dN}{dE} = \frac{1}{p_\pi},\qquad E_\pm = \frac{1}{2}(E_\pi \pm p_\pi)$$

The decay asymmetry $\alpha = (E_1 - E_2)/(E_1 + E_2)$ is uniformly distributed as $dN/d\alpha = \beta/2$ with limits $\alpha_\pm = \pm\beta$. An asymmetry cut of $\alpha_{max} = 0.8$ retains 80% of decays while drastically reducing the combinatorial background.

### $\pi^0$ Invariant Mass Measurement

Photon pairs with an invariant mass near the $\pi^0$ mass ($m_{\pi^0} = 138 \pm 10$ MeV, 3$\sigma$ window) are identified as $\pi^0$ candidates. The $\pi^0$ mass peak is fitted with a Gaussian for the signal and a second-order polynomial for the background, with results produced separately for each sector of the EMCal.

### Associated Charged Hadron Acceptance

A bootstrap method calculates the hadron acceptance correction using the charged hadron yield previously measured by PHENIX (PPG023). The published yield is fitted with a cubic spline at low $p_T$ and a power law $a p_T^{-k}$ at high $p_T$, connected by a Woods-Saxon at $p_T = 3$ GeV. The expected hadron level in each $\xi$ and $p_{Tt}$ bin is calculated by integrating the fit over the appropriate $p_{Th}$ range weighted by the trigger $p_T$ distribution. A centrality-independent systematic error of 7.5% from the published hadron measurement is attributed to the efficiency.

### Background Term Extraction

Background terms are extracted from the acceptance-corrected mixed event per-trigger yields:

$$n_t Y_{mix}^{\psi} = 4\pi \langle a_t \rangle \langle a_h \rangle + 2\pi \langle \sigma b_t \rangle \langle \sigma b_h \rangle \cos(2\Delta\phi)$$

$$n_t^{\pm} Y_{mix}^{\psi\pm} = 2\pi \langle a_t \rangle \langle a_h \rangle \pm 4 \langle a_h \rangle \langle b_t \sigma \rangle \pm 4 \langle a_t \rangle \langle b_h \sigma \rangle \cos(2\Delta\phi) + \pi \langle b_t \sigma \rangle \langle b_h \sigma \rangle \cos(2\Delta\phi)$$

The integral of $n_t Y_{mix}^{\psi}$ reveals the flat term while a $\cos(2\Delta\phi)$ fit measures the cosine amplitude. The in/out-of-plane difference $n_t^+ Y_{mix}^{\psi+} - n_t^- Y_{mix}^{\psi-}$ isolates the cross terms.

### Background Component Isolation

Components are extracted for use in $\kappa$ corrections:
- $\langle a_t \rangle$: counting triggers in bins of $p_{Tt}$ and centrality
- $\langle a_h \rangle$: dividing the flat background term by $\langle a_t \rangle$
- $\langle b_h \sigma \rangle$: from the cosine amplitude specific to the in/out-of-plane analysis
- $\langle b_t \sigma \rangle$: from the additional flat background in in/out-of-plane yields
- $\langle b \rangle \approx \langle b \sigma \rangle / \langle \sigma \rangle$: used for estimating $v_2$

### Background Subtraction

The $\kappa$-corrected backgrounds are subtracted from the same-event per-trigger yields in a stepwise procedure. For the all-plane case, three steps are shown: raw data with flat background overlaid, intermediate result after flat subtraction with the cosine background overlaid, and final result with both backgrounds removed. For the in/out-of-plane case, five steps are required due to the additional cross terms. The $4 \kappa_{ht} \langle a_h \rangle \langle b_t \sigma \rangle$ flat term and $4 \kappa_{th} \langle a_t \rangle \langle b_h \sigma \rangle \cos(2\Delta\phi)$ term flip sign between in and out-of-plane, as required by the derivation.

### $v_3$ Contribution

A possible $\cos(3\Delta\phi)$ modulation is observed in the per-trigger yields, attributed to triangular flow ($v_3$). The single particle distribution expanded to arbitrary flow harmonics:

$$\frac{dn}{d\phi} = a + \sum_{n=2}^{\infty} b_n \cos\big(n(\phi-\psi_n)\big)$$

The pair distribution becomes:

$$n_{th} = 2\pi a_t a_h \left(1 + 2 \sum_{n=2}^{\infty} v_{tn} v_{hn} \cos(n\Delta\phi)\right)$$

Using $v_3$ estimates from previously published PHENIX results, the $\cos(3\Delta\phi)$ contribution is well described. By integrating the away-side yield over $2\pi/3 < \Delta\phi < \pi$, the $v_3$ contribution is completely removed. This motivates the choice of integration limits throughout the analysis.

### Results

Per-trigger $\pi^0$ - $h$ yields from 0--60% Au+Au at $\sqrt{s_{NN}} = 200$ GeV are produced for four trigger $p_T$ bins and three $\xi$ bins, consistent with PHENIX preliminary results. The ratio of central (0--20%) to mid-central (20--60%) integrated away-side yields reveals a significant suppression, indicating path-length-dependent energy loss.

### Reaction Plane Dependence and Away-Side $v_2$

Assuming the azimuthal distribution of associated particles:

$$j(\phi_a) = p(\phi_a - \phi_t) \Big(1 + 2 v_{2a} \cos\big(2(\phi_a - \psi_2)\big)\Big)$$

and taking the ratio of in/out-of-plane to all-plane yields before integration removes the unknown associated particle distribution $p(\Delta\phi)$:

$$\frac{Y^+(\Delta\phi) - Y^-(\Delta\phi)}{Y(\Delta\phi)} = \frac{16\sigma\, v_{2a} \cos(2\Delta\phi) (a_t^2 - b_t^2)}{\pi (a_t^2 - 4\sigma^2 b_t^2/\pi^2)(1 + 4 v_{2a} b_t \cos(2\Delta\phi) / a_t)}$$

For small trigger $v_2$:

$$\frac{\pi}{16\sigma} \frac{Y^+(\Delta\phi) - Y^-(\Delta\phi)}{Y(\Delta\phi)} \approx v_{2a} \cos(2\Delta\phi)$$

For integrated yields rather than full $\Delta\phi$ distributions:

$$\frac{\pi}{16\sigma} \frac{\int_{2\pi/3}^{\pi} Y^+(\Delta\phi) - \int_{2\pi/3}^{\pi} Y^-(\Delta\phi)}{(0.75 \pm 0.1) \int_{2\pi/3}^{\pi} Y(\Delta\phi)} = v_{2a}$$

A positive away-side $v_2$ is measured for all $\xi$ bins, decreasing at higher trigger $p_T$.

### Systematic Errors

Three sources of systematic uncertainty are quantified:
1. **Hadron efficiency**: ~7.5%, centrality-independent, from the published hadron yield (does not affect ratios)
2. **Background subtraction**: including $\kappa$ correction uncertainties, propagated through the subtraction procedure
3. **Momentum resolution**: studied via fast Monte Carlo simulations ($\sigma_p/p \simeq 1\% \oplus 1\%p$), with the relative difference between smeared and unsmeared integrated away-side yields taken as the systematic error

---

## Direct Photon Tagging

### Separation of Direct and Decay Photons

To separate direct photons from the inclusive sample, a tagging method is implemented. Photons from particle decays are tagged as *decay* while all others (fragmentation, bremsstrahlung, and prompt) are tagged as *direct*. The tagged photon counts form a linear system:

$$\begin{pmatrix} N_{\text{dec}}^t \\ N_{\text{dir}}^t \end{pmatrix} = \begin{pmatrix} \tau & \mu \\ 1-\mu & 1-\tau \end{pmatrix} \begin{pmatrix} N_{\text{dec}} \\ N_{\text{dir}} \end{pmatrix}$$

where $\tau$ is the tagging efficiency and $\mu$ is the mistag fraction. Inverting:

$$N_{\text{dir}} = \frac{\tau}{\tau - \mu} N_{\text{dir}}^t - \frac{1-\tau}{\tau - \mu} N_{\text{dec}}^t$$

### Mistag Probability

The mistag probability (tagging a random photon as decay) is obtained via simulations. Photons are simulated flat in $p_T$, $\phi$, and $\eta$ and embedded into real PHENIX event data. Reconstructed clusters passing the photon selection cuts are passed through the decay tagging algorithm. The ratio of decay-tagged candidates to all candidates gives $\mu(p_T, \text{centrality})$. The distributions are fit with functions of the form $p_0/p_T + p_1/p_T^2 + p_2/p_T^3$.

### Tagging Efficiency

$\pi^0$'s are simulated flat in $p_T$, $\phi$, $\eta$ and forced to decay to $2\gamma$. 80% of decay photons are produced via this channel; the remaining 20% are primarily from $\eta \rightarrow 2\gamma$. Both are simulated, embedded, and run through the same reconstruction and tagging procedure. The composite tagging efficiency:

$$\tau(p_T) = \frac{0.8 \left(N_{\text{tag}}^{\pi^0}/N_{\text{true}}^{\pi^0}\right) + 0.2 \left(N_{\text{tag}}^{\eta}/N_{\text{true}}^{\eta}\right)}{0.8 \left(N_{\text{meas}}^{\pi^0}/N_{\text{true}}^{\pi^0}\right) + 0.2 \left(N_{\text{meas}}^{\eta}/N_{\text{true}}^{\eta}\right)}$$

The efficiency peaks at $\sim 10$ GeV where the shrinking partner search window increases acceptance, and falls at higher $p_T$ due to cluster merging effects.

### $R_\gamma$ Measurement

The ratio of inclusive to decay photons:

$$R_\gamma = 1 + \frac{N_{\text{dir}}}{N_{\text{dec}}}$$

is measured using the tagging method and is consistent with previously published PHENIX results for all centralities. The in vs. out-of-plane $R_\gamma$ ratio is most often less than one, consistent with a positive decay photon $v_2$.

---

## Direct Photon - Hadron Correlations

### Initial Tagging Method

The tagged photon - hadron yields can be written as a linear combination of the true yields:

$$\begin{aligned}
N_{\text{dec}}^t Y_{\text{dec}}^t &= \tau N_{\text{dec}} Y_{\text{dec}} + \mu N_{\text{dir}} Y_{\text{dir}} \\
N_{\text{dir}}^t Y_{\text{dir}}^t &= (1-\mu) N_{\text{dir}} Y_{\text{dir}} + (1-\tau) N_{\text{dec}} Y_{\text{dec}}
\end{aligned}$$

The per direct photon trigger yield:

$$Y_{\text{dir}} = \frac{\big(\frac{\tau}{\tau-\mu}\big) N_{\text{dir}}^t Y_{\text{dir}}^t - \big(\frac{1-\tau}{\tau-\mu}\big) N_{\text{dec}}^t Y_{\text{dec}}^t}{\big(\frac{\tau}{\tau-\mu}\big) N_{\text{dir}}^t - \big(\frac{1-\tau}{\tau-\mu}\big) N_{\text{dec}}^t}$$

The $\xi$ binning renders the parent $p_T$ bias of the tagging method negligible (< 10% effect, covered by statistical errors). While the measurement is successful, the lack of precision caused by the relatively low statistics of the tagged decay photon sample (tagged direct photons often exceed tagged decay photons by a factor of 10) limits the significance of the results.

### Hybrid Method

To overcome the statistical limitations of the pure tagging method, a hybrid approach is developed. The per decay $\gamma$ - $h$ integrated away-side yields are obtained from the per $\pi^0$ - $h$ yields by sampling the measured $\pi^0$ trigger integrated away-side yield as a function of $p_{T,\pi^0}$ while weighting by the parent $p_T$ distributions of simulated decay photons. The per inclusive $\gamma$ - $h$ yields undergo the same background subtraction procedure as the $\pi^0$ - $h$ yields. The direct photon - hadron yield is then:

$$Y_{\text{dir}} = \frac{1}{N_{\text{dir}}} \big(N_{\text{inc}} Y_{\text{inc}} - N_{\text{dec}} Y_{\text{dec}}\big)$$

where $N_{\text{dir}}$ and $N_{\text{dec}}$ are obtained from the tagging method ($R_\gamma$).

### Results

The direct photon triggered integrated away-side yield in 0--60% collisions is consistent with PHENIX PPG210. The ratio of central to mid-central yields reveals a several-sigma suppression, approximately constant in $\xi$ and $p_T$. The away-side $v_2$ of hadrons associated with direct photon triggers is consistent with zero, with improved statistical precision over the pure tagging method.

---

## Discussion

The method requires no external $v_2$ measurements — all flow-dependent backgrounds are handled by the reaction-plane-binned event mixing. $v_3$ is eliminated via integration limits. Results are consistent with published PHENIX measurements, validating the approach. Central-to-mid-central suppression and positive away-side $v_2$ point to path-length-dependent energy loss in the QGP. The suppression is similar for $\pi^0$ and direct photon triggers despite differing surface bias and parent parton composition.

---

## Code Structure

```
src/                          C++ analysis framework
  gamAna.h/cc                 SubsysReco module: event loop, track/photon handling
  gamAnaTree/                 Per-event tree building and histogram filling
  EmcCluster.h/cc             EMCal cluster reconstruction and shower shape cuts
  CentralTrackData.h/cc       CNT track matching, DCA, and PID (e/p, mass^2)
  TriggerParticle.h/cc        Trigger selection (photon, pi0, asymmetric cut)

macros/                       ROOT/C++ analysis and plotting (~70 macros)
  artist_gamAnaTree.C         Main driver: builds dphi histograms, applies mixing
  acceptCorrect.C             Acceptance function from RP-dependent mixing
  calcKComponents.C           kappa corrections
  subtractKWeightedBG.C       Full background subtraction
  integAwaySide.C             Integrated yield extraction and v2
  CalcRPReso.C                Reaction plane resolution from sub-events
  getFlow.C                   Flow background isolation
  draw*.C                     Publication-quality figures
```

---

Z. Rowan, "Removal of Anisotropic Background from Neutral Pion and Tagged Direct Photon–Hadron Correlations of Au+Au 200 GeV Collisions," PhD thesis, CUNY Graduate Center, 2022.
