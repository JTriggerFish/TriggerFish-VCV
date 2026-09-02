# Electric Piano modelling notes

This is a working engineering record for the sample-free Electric Piano in
`src/models/ElectricPiano.hpp`. It records the provenance and limits of the
current real-time model so that listening calibration cannot quietly turn into
claims of circuit or mechanical accuracy. It is a technical-report stub, not a
finished validation report.

## Evidence labels

- **Source**: directly supported by a cited service document, measurement, or
  published model.
- **Derived**: follows mathematically from a sourced topology or component
  value, subject to stated assumptions.
- **Calibrated**: selected by rendered-signal tests and listening within the
  model's physical structure.
- **Provisional**: plausible but awaiting measurements from an identified
  Rhodes assembly.

## Signal path and rates

The note path is hammer/action -> tine and tone bar -> magnetic pickup. Sixteen
voices and released-note tails are summed into one preamplifier -> tone network
-> optical stereo router -> two power modules -> electrical speaker loads. The
electrical load affects output-device current and shared supply sag, but its
cone motion is not rendered as an acoustic cabinet simulation.

| Section | Rate | Status |
| --- | ---: | --- |
| Coupled tine/tone-bar modes | host | reduced physical model |
| Hammer contact while engaged | 64x host | nonlinear physical model |
| Magnetic pickup | 4x host | precomputed calibrated 2D flux field |
| Figure 11-8 input transistor pair | 2x host | five-unknown nonlinear circuit solve |
| Figure 11-8 active tone feedback | 2x host | thirteen-unknown nonlinear circuit solve |
| Figure 11-8 volume/optical routing | 2x host | linear audio-path reduction |
| Figure 11-9 power modules | 2x host | coupled transistor-level circuit solve |
| Control smoothing and shared supply recovery | host | calibrated slow states |

Seventh-order Chebyshev interpolation/decimation surrounds the 4x nonlinear
pickup path and a separate 2x amplifier path. The contact integrator substeps
only while contact exists. Modes taper before host Nyquist instead of being
allowed to alias.

## Action and velocity

`Vel Curve` applies a controller-response curve. `Dynamics` maps the result
through one explicit action-law calibration (`velocity^0.85`) to the published
0--4 m/s incoming hammer-speed range. A very small positive velocity is
treated as controller residue before a collision starts; an active collision
is never terminated by a timer. **Calibrated:** the default Dynamics value is
the full input range. The present action is not a multi-lever key escapement
simulation, so hammer speed as a function of key travel remains an open model.

## Hammer contact

The finite hammer mass and all participating modal velocities are coupled in
both directions. For compression `q > 0`, the contact force has the
Hunt-Crossley form

`F = k q^n (1 + alpha qDot)`, clamped to a non-tensile force.

The solve uses the published quantities directly: 11 g hammer mass, 4 m/s
maximum speed, `n = 2.8`, `k = 1.5e11 N/m^n`, and damping weight
`lambda = 9e10 N s/m^(5/2)`. `lambda/k = 0.6 s/m` is retained when Hammer
moves stiffness around the sourced middle-C anchor. Hammer mass and exponent
no longer change with a panel knob. The service manual's five tip zones (30,
50, 70, 90 durometer, then wrapped extra-hard) select the factory baseline;
the panel span and mapping from durometer to stiffness remain clearly named
listening calibrations. Separation follows zero compression or the non-tensile
limit of the dissipative law, never a timer.

At default middle C, rendered contact duration changes only modestly with
velocity. The nonlinear pickup and contact force together provide the larger
low-harmonic change; the first inharmonic attack coordinate rises more gently
and is no longer used as a proxy for all velocity brightness.

The retained seven attack coordinates are projected over a five-point spatial
integration of the finite neoprene contact footprint. Falaize's downloadable A4
model uses a 15 mm distributed-force width; the real-time calibration uses a
conservative 6--12 mm effective loaded strip from hard to soft. A trial deriving
this longitudinal dimension from Hunt--Crossley indentation was removed: it
assumed an unsupported spherical 4 mm tip even though the Rhodes tips are
block-shaped, and narrowed the contact enough to disturb the validated bass.
Velocity affects temporal contact bandwidth through the nonlinear force law,
not through an invented velocity-dependent tip geometry. There is no
velocity-cubed modal multiplier or post-contact brightness envelope.

Pickup bandwidth and contact bandwidth are deliberately separate. All seven
defined attack coordinates continue to load the hammer in the 64x collision
solve even when an upper key places a coordinate above the 48 kHz pickup
model's Nyquist limit. Such a coordinate is excluded from pickup reconstruction
and discarded only after physical separation. The former shared activity flag
removed most higher-mode point mobility from treble collisions, concentrating
the force spectrum in the last audible coordinate. The corrected split leaves
the bass bit-for-bit unchanged where all seven coordinates were already in
band, while redistributing the E5 onset across the first two audible attack
coordinates.

An offline uniform-beam tail check does not justify adding more coordinates
yet. With the present finite footprint, the next eight ideal modes contribute
about 2--7% as much unweighted impulse mobility as the retained seven at E5--E6,
before temporal contact bandwidth attenuates them further. In the bass that
upper bound is larger, but the current sound and regression anchors are already
successful. Extending the tail would first require a converged continuous
pressure-footprint integral: the current five-point calibrated pressure profile
can spatially alias modes beyond those it was designed to approximate.

Trials that directly boosted upper modal contact projection were rejected.
Because the nonlinear hammer sees the sum of modal point mobilities, changing a
projection also changes contact impedance and moves zeros in the hammer-force
spectrum; the same boost suppressed E5 while exaggerating E6. This is evidence
of an underconstrained upper-tip/contact calibration, not evidence against the
Hunt--Crossley topology. Force histories and pressure footprints for the five
tip grades remain the measurements needed to close that gap.

**Source:** Hunt and Crossley's nonlinear damping term, the Rhodes hammer-tine
topology, the numerical fit and material data reported by Sonderbo, and the
service-manual tip zones. **Derived:** bidirectional modal force, SI modal
coordinates, and finite-width mode-shape projection. **Calibrated:** MIDI-to-
hammer action law, durometer-to-stiffness mapping, panel stiffness span, patch
width, five-point pressure profile and keyboard scaling. Force/displacement or
contact-duration measurements from tips of known hardness and a measured tip
profile are still needed for absolute validation.

## Tine, tone bar and mounting block

The fundamental is a two-coordinate asymmetric fork. Tine and tone-bar mass,
stiffness and common-base coupling form a 2x2 generalized eigenproblem. The
two prongs' uncoupled fundamentals are 0.997--0.998 apart, consistent with the
service manual and SLDV observation that each prong carries a strong played
fundamental. Coupling splits them into two near-unison normal coordinates; the
played coordinate is calibrated back to the requested pitch. Normal
frequencies, modal masses, strike projections and tine pickup projections all
follow from the same eigenvectors.

The tone bar's **sub-fundamental is a separate mode**, not the lower member of
that fundamental pair. This corrects an earlier structural mistake: detuning
the tone-bar fundamental to 0.45--0.75 of the tine manufactured a persistent
submode and pickup-sideband family whose register trend was opposite the SLDV
data. Gabrielli/Cantarini measure the separate submode at 0.83*f0 and -9.1 dB/s
for F1, and 0.58*f0 and -138 dB/s for F3. Its initial displacement residues are
-28 and -4.5 dB relative to the played mode. The model passes through all six
anchors. One tine-participation coefficient projects both hammer force and
pickup observation, so the fit remains reciprocal rather than adding an
output-only oscillator. Above F3, frequency, decay and participation follow a
smooth provisional fit to aligned direct-harp recordings; by the upper register
the submode is about 40--60 dB below the played mode during the first 100 ms.
The magnetic pickup, not the mechanics, generates its measured sidebands.

The separate submode shares the fork's common-base transfer. `Coupling` scales
its tine participation from the played fundamental pair's tone-bar eigenvector,
with a small distributed mounting path, and applies the same scale at hammer
projection and pickup observation. The exact panel midpoint remains the
published-residue calibration; the two endpoints change submode residue by
roughly -6/+3 dB rather than leaving Coupling as an almost exclusively late-
decay control.

Mounting loss follows the squared common-base bending reaction, normalized by
modal mass. It approaches zero for the played anti-phase mode as tine and
tone-bar forces cancel. Coupling therefore changes hybridization and decay
rather than crossfading two oscillators.

The higher-mode model retains seven signed attack coordinates. Gabrielli and
Cantarini's measured F1 ratios 7.2/20.6 and F3 ratios 7.4/20.7/38.7 are now
explicit keyboard anchors. A Rayleigh--Ritz cantilever with a sliding point mass
fits those anchors with tuning-mass/beam-mass fractions 0.146 and 0.177 and
positions 0.859L and 0.840L. The first fit uses the two published F1 modes; the
F3 fit minimizes the residual over all three published modes. Extrapolating the
mass fraction versus published tine length and the spring position versus key
produces the reduced treble curves implemented in the real-time model. Their
top-key pre-shear ratios are approximately 7.791/20.599/38.144, rather than the
former unloaded-beam 6.267/17.548/34.387 values. That former endpoint came from
Falaize and Helie's A4 Euler--Bernoulli example, not a measurement of a
production A4 assembly; every real Rhodes tine still carries its tuning spring.

The point-mass eigenproblem is fitted offline and reduced to smooth ratio
polynomials so continuous pitch modulation does not solve a generalized matrix
problem per sample. Earlier versions stopped there: point-mass-loaded
frequencies were combined with unloaded cantilever shapes and the uniform
quarter-beam modal mass. That is not one eigenproblem and gave internally
inconsistent strike and pickup projections.

Modes 3--5 now use the corresponding spring-loaded Rayleigh--Ritz eigenvectors
and generalized masses. Each tip-normalized loaded shape is expanded in the
first eight uniform cantilever shapes; its eight coefficients and modal-mass
ratio are reduced to fifth-order keyboard polynomials. Against the offline
73-key solve, maximum coefficient error is below `6.5e-6` and mass-ratio error
below `2.3e-5`. Contact-footprint integration and pickup observation therefore
use the same loaded shape as the frequency fit, while the real-time cost remains
a handful of polynomial evaluations at note-on/control refresh.

This reduction is **derived/provisional** above F3 because the paper's claimed
whole-key modal dataset is not present in the openly accessible material. The
already tuned bottom-key 7.11/20.25 pair is retained over the one-semitone
approach to F1, and F1/F3 remain exact. Only the three modes constrained by
measurements use the point-mass curves. Applying unconstrained high-order
point-mass eigenvalues would alter the good lower-middle attack, while those
modes rapidly approach the audible-band taper; they retain the distributed-beam
ratios pending measured residues.

Beyond the final F3 measurement, a normalized first-order
Timoshenko/Rayleigh correction additionally accounts for shear and rotary
inertia in progressively shorter tines. It leaves the fundamental and all
measured bass anchors exact, but lowers the first three attack-mode ratios at
the top key by approximately 0.44%, 1.36% and 2.69%. This removes 8--47 cents
of ideal-beam error without inventing another resonant family.

The same F1/F3 fits report both modal magnitude and amplitude-decay slope:
21.1/67.7 dB/s at F1 and 294/37/161 dB/s at F3. The slopes convert to
exponential lifetimes with `tau = 8.6858896 / |slope|`. Mode-specific smooth
multipliers make Sonderbo's `sigma0 + sigma1 k^2` distributed-loss law pass
through those measured values. F3 is the final measured damping anchor, so its
multiplier is not extrapolated unchanged through the treble. Above F3 it relaxes
independently by A4 and reaches the sourced distributed-loss law. It is not tied
to the point-mass frequency curve: a tuning spring can shift a mode without
preserving F3's anomalous loss. This
corrects a genuine extrapolation error which left the fourth and fifth attack
modes about eleven and five times too long in the upper register. At A4 their
amplitude lifetimes now follow the distributed-loss law at roughly 8.9 and
4.5 ms instead of 99 and 23 ms. The first attack mode simultaneously falls from
the F3 fit back to its 25 ms beam lifetime, retaining a quieter short-lived
presence while the upper modes supply a briefer onset sparkle.

Table 8.6 also reports the displacement magnitude of the first inharmonic mode
as -65 dB against a -10 dB F1 fundamental and -59.5 dB against a -4.5 dB F3
fundamental: -55 dB relative at both anchors. A trial observation-only factor
forced the analytic coordinate through those values, but it broke reciprocity,
removed the previously successful bass attack and did not remove the offending
sideband because that came from the incorrectly constructed fork mode. It has
been removed. The spring-loaded eigenfunction and generalized mass remain in
both force and observation paths. Absolute 7x residue is again an explicit
calibration gap: matching an SLDV displacement line and matching a pickup
recording cannot be conflated, because the nonlinear pickup also puts a strong
integer/near-integer family in the same spectral region.

The same table's 20x and 39x magnitudes are recorded as calibration targets but
are not forced independently yet. In the present finite-duration contact solve,
the F3 39x coordinate lies close to a force-spectrum null; matching it with an
output residue alone requires roughly 483x gain. Testing that shortcut broke
reciprocal gain staging and alias headroom and regressed the bass. Those two
residues need a joint fit of measured hammer-force/contact bandwidth and modal
observation rather than a large hidden multiplier.

The relaxation coordinate is **derived/provisional**, not another measurement:
published damping values stop at F3. Freezing the F3 anomaly indefinitely has
less physical support than returning to the sourced distributed beam-loss law,
but treble decay measurements should eventually replace that assumption.

Applying the decay columns alone was a calibration error: it made the F3 3.6
and 6.8 kHz modes roughly eleven and five times longer without accounting for
the fitted residues, greatly increasing their total energy. Generalized
contact projections are therefore normalized by `sqrt(tau_beam / tau_fit)` so
the measured lifetime update preserves the previous impulse-energy calibration.
This is a modal normalization used by the single model at every key, not a
crossfade between old and new instruments. Rendered modes use a fixed 48 kHz
pickup bandwidth at host rates of 48 kHz and above, so higher host rates cannot
introduce different audible modal content. Contact-only modes instead use the
64x collision bandwidth and never enter pickup output.
The former 15.36--21.6 kHz taper ran before the 4x pickup and discarded valid
short attack modes unnecessarily. Full gain now extends to 19.2 kHz and tapers
to zero at 23.52 kHz; the seventh-order 4x decimator remains responsible for
the final anti-alias filtering.

`Bell = 0.52` is unity gain for the calibrated attack-mode residues. Full panel
travel spans 24 dB in the bass/middle and grows smoothly to 36 dB by C6, where
fewer mechanical coordinates remain below pickup bandwidth. The same gain
tracks their signed transverse pickup weights. Above middle C, Bell additionally
moves the pickup's vertical observation point over a 0.09 mm full span, blended
smoothly to full travel at C6. This is the same physical voicing dimension by
which a technician changes harmonic bell at the pickup, rather than enormous
gain on an almost absent treble attack coordinate. The combined Bell/Tone
position is bounded by Tone's existing physical alignment range. Trajectory
normalization compares the shifted position with the same key's neutral Bell
trajectory, so even extreme Tone and Proximity combinations remain bounded.
This is a sound-design balance around the physical default, not a hidden default
correction or a broadband treble shelf. The 80-note default render corpus is
bit-identical because the added offset is exactly zero at 0.52.

`Decay = 0.5` is exactly a unity lifetime multiplier, not a compensating
envelope setting. At default Coupling, the played upper fork mode has amplitude
time constants of about 3.7 s at the bottom of the 73-key range, 1.2 s at C4,
and 0.9 s at the top (approximately 26, 8.3 and 6.2 s respectively to -60 dB
for a pure exponential). This decreasing sustain toward the treble is retained
as the neutral calibration. The higher-mode losses are sourced as described
above; the fundamental mounting-loss scale remains a provisional calibration
until long free-decay measurements across a serviced harp are available.

Pfeifle includes shear-beam and large-deflection coupling between two
polarizations. The shear correction is now retained as described above, and the
real-time reduction already has a weak transverse body coordinate and a
two-dimensional pickup trajectory. A cubic large-deflection force was assessed
again after restoring the upper pickup harmonic bed, but is still not added:
no published coefficient or keyed tine trajectory is available to calibrate
that term. Adding one now would spend CPU and risk recreating a dissonant attack
without addressing a measured deficit.

The Strike midpoint follows a curved keyed line rather than the former straight
0.40L--0.20L approximation. Its explicit checkpoints are 0.38L at the bass
end, 0.29L at C3, 0.22L at F3, 0.205L at C4, 0.16L at C6 and 0.14L at the top.
Smooth interpolation between them removes the F3-region fundamental
cancellation and excessive middle/treble attack-mode projection found in the
73-key render sweep. This shape mirrors the service manual's instruction to
set C4, F3 and C3 separately at the clearest, maximum-power strike point, then
accept the intervening keys. The dimensions remain calibrated rather than
claimed as measured factory jig settings. The panel now moves a signed physical
striking-line offset around that midpoint: full travel tapers smoothly from
plus/minus 6 mm in the bass to plus/minus 1 mm at the top. It no longer maps the
two halves independently toward 0.04L and 0.96L, whose centre slope differed by
up to 8.2x and crossed upper-mode nodes abruptly. A 73-key regression checks
equal left/right local derivatives, bounded physical travel and inter-key
continuity. **Source:** the Rhodes service manual's unequal-leg
tuning fork and striking-line procedure, Sonderbo's register direction, and
the separated-prong/modal measurements cited below. **Derived:** the
coupled-fundamental eigenproblem, mass orthogonality and bending-support
reaction; reciprocal projection of the separate tone-bar submode.
**Calibrated/provisional:** keyed strike checkpoints, millimetre panel span,
upper-key submode continuation, tuning-spring perturbations between measured
anchors, contact width and fundamental loss scaling. Laser-vibrometer and
isolated-prong data across several keys would supersede these calibrations.

## Magnetic pickup

Two transverse tine coordinates move through an asymmetric magnetic field.
Pfeifle derives a three-part integral over an idealised pole cross-section,
but their paper supplies a schematic rather than production dimensions and
explicitly lists geometrically accurate pickup modelling as future work. A trial
implementation guessed those dimensions and produced an excessive second
harmonic, so it has been removed rather than hidden behind a register blend.
Falaize and Helie's published finite circular-pole law was also evaluated as a
whole-keyboard replacement. It is useful in their measured A4 comparison, but
with this model's keyed physical excursions it raised lower-register H2 by
roughly 5--14 times, about tripled peak level, and further depleted H4+ in the
treble. That fails the lower-register reference and demonstrates why an A4
pickup reduction cannot be presented as measured Rhodes pole geometry.
The model uses the previous listening-calibrated pair of softened flux edges,
clearly classified as a reduced field rather than measured geometry. A 241 by
135 table stores its two-dimensional vertical/radial gradient; bilinear lookup
at the 4x pickup rate is checked against direct evaluation to 0.3% relative
error. Horizontal tine motion changes radial clearance by the chain rule, and
voltage remains the time derivative of one scalar flux construction.

`Tone` moves the vertical resting alignment; the upper-register part of `Bell`
makes the smaller bounded movement described above; `Proximity` changes the
pole gap.
The gap is in millimetres: the panel covers the manual's 1/8-inch ordinary
clearance down to the documented 0.020-inch close setting. The neutral per-key
screw position stays exactly 1.6 mm through middle C, preserving the calibrated
bass and bark, then moves smoothly to 0.52 mm at the top key, just above the
manual's close setting. The former 0.60 mm endpoint left too little pickup-
curvature contribution in the highest register. This is a physical pickup gap,
not a register EQ. Proximity moves
logarithmically around that keyed neutral position. The closest endpoint remains
zoned so bass keys retain the ordinary 1/16-inch limit.

Small-signal alignment normalization is supplemented by the geometric-mean RMS
of three representative elliptical tine trajectories. This one scalar per
pickup setup preserves instantaneous curvature and overload while preventing
Proximity from becoming another level control. At velocity 0.86,
default-to-closest renders across seven register points stay within about
-0.81 to +0.48 dB, instead of cutting bass while boosting treble by roughly 6 dB.

The keyed neutral gap is a spectral calibration rather than a level effect.
Before it, default H4-and-above attack energy was about -36 dB at C6 and -61 dB
at C7 relative to the fundamental, leaving valid mechanical lines perceptually
isolated. It is now about -22 dB and -38 dB respectively. The complete harmonic
family exceeds the measured attack-mode bins in the guarded C6, C7 and top-key
renders. The pickup-gap path is unchanged through middle C; the captured
F1, E2, C3 and F3 spectral, dynamic and envelope fingerprints remain inside
their regression bounds after the independent high-audio modal-taper update.

Nineteen explicit, four-key pickup-level checkpoints stand in for the
individual screw adjustment performed on a real harp; a bass-only excursion
correction documents where the reduced uniform-beam model still over-predicts
long-tine travel. The checkpoint fit uses the geometric-mean RMS of
velocities 0.30, 0.60 and 0.90 and preserves the middle-C amplifier voltage
calibration. Over all 73 keys it has a deliberate -0.42 dB/octave voltage tilt,
a 2.46 dB end-to-end span, and no adjacent-semitone step above 0.25 dB. The
post-C4 checkpoints were re-fitted after the point-mass update; all lower
checkpoints remain unchanged.

The seventh-order 4x decimator is now the only pickup reconstruction filter.
The removed host-rate 16.5 kHz one-pole was not present in the pickup or
Peterson schematic and attenuated about 0.6 dB at 8 kHz, 1.1 dB at 12 kHz and
1.9 dB at 20 kHz. Register-wide folded-spectrum tests remain responsible for
alias rejection instead of using that unsourced pole as a safety blanket.

**Source:** published Rhodes field topology and the service adjustment
procedure. **Derived:** scalar-flux differentiation and velocity-dependent
harmonic/intermodulation generation. **Calibrated:** reduced edge-field shape,
excursion scale, per-key gap, trajectory compensation and sensitivity. These
trims are kept separate from published mechanical constants. A measured
two-dimensional flux map and pickup impedance are still needed. The model does
not contain a separate bark waveshaper.

## Damper and mechanical noise

The pitched damper is a modal loss change, separate from the audible damper
mechanics event. Sonderbo's published damped-spring connection uses
`K1 = 100 N/m`, `K3 = 100000 N/m^3`, and `R = 0.5 kg/s`. Its strongly
overdamped slow relaxation is approximately `R/K1 = 5 ms`. The reduced modal
model uses that value as the fast endpoint of a logarithmic 5 ms--1.2 s Release
range. The default `Release = 0.12` gives a 9.65 ms amplitude time constant
(about 67 ms to -60 dB); rendered key-up output is already about 55 dB down in
the 50--100 ms window across the keyboard. This is consistent with the service
manual's requirement that a properly tensioned felt damp the sound
"immediately", while leaving the rest of the knob for deliberate long-release
design. The paper's values are model parameters rather than measurements of a
particular restored instrument, so `ElectricPianoDefaultRelease` remains a
single visible listening checkpoint.

The current reduction does not solve the published cubic damper connection as
a second nonlinear contact. It applies its calibrated relaxation to each modal
coordinate; isolated release recordings would be required to justify that
extra state and cost.

Key-down, key-release and damper contacts excite independent short mechanical
resonators. Every event belongs to its own voice/tail and begins at that MIDI
event, avoiding the earlier shared random burst that lost timing on chords.
Mechanics is mixed outside the magnetic pickup but before the shared amplifier.
The default is deliberately subordinate to pitched energy. **Calibrated:** all
noise spectra, levels and durations; isolated keybed and damper recordings are
needed for measurement-based fitting.

## Peterson preamplifier and controls

Figure 11-8 and Figure 11-9 must not be conflated.
`src/models/PetersonPreAmplifier.hpp` now solves the audio path's nonlinear
input section as a five-unknown transistor circuit:

- Q1 is the selected 2N3392 common-emitter input device, with the drawn 390 kOhm
  rail bias, 33 kOhm ground bias, 12 kOhm collector and 1.5 kOhm emitter paths;
- Q2 is the directly coupled selected-2N3392 emitter follower with its 4.7 kOhm
  emitter path;
- the 0.22 uF input capacitor, 330 pF Q1 base/collector capacitor and 5 uF
  output capacitor are trapezoidal companions in the same Newton solve;
- the nominal 64 kOhm following load represents the parallel loading visible at
  the entrance to the tone network.

This is a circuit reduction, not a preamp waveshaper. Both devices use Ebers-
Moll terminal currents and an analytic Jacobian. At 250 Hz, 0.2 V peak, the
checked-in ngspice reference gives 1.06733 V RMS and 0.07832% THD; real time
gives about 1.0664 V RMS and 0.07815%. At 1.5 V peak, ngspice gives 7.68424 V
RMS and 5.2267% THD; real time gives about 7.68 V RMS and 5.1902%. The strong
overload is asymmetrical and dominated by low harmonics because it follows from
Q1 cutoff/saturation and Q2 headroom, not a fitted clip curve.

`src/models/PetersonTonePreAmplifier.hpp` now solves Figure 11-8's active tone
section as the circuit drawn, rather than as the former pair of one-pole shelves.
The 100 kOhm Bass and Treble controls are each represented by their two moving
pot segments. Their 0.0047, 0.047 and two 0.1 uF branches, 22 kOhm bridge,
6.8 kOhm end arms and 390 kOhm/1 MOhm feedback divider all participate in the
same zero-delay solve as Q3/Q4 and Q5. Q3/Q4 are the rotated,
Darlington-connected selected-2N3392 common-emitter stage visible in the service
drawing; Q5 is the following compensated selected-2N3392 emitter follower.
Consequently the tone positions change both the small-signal response and the
stage's nonlinear feedback/headroom instead of feeding an unrelated waveshaper.

The checked-in ngspice deck now covers the input pair and this complete active
tone section. At 250 Hz it produces 4.4299 mV RMS for 1 mV peak input, 2.2051 V
RMS for 0.5 V peak and 4.2289 V RMS with 5.34% THD for 1 V peak. The split
real-time solves produce 4.4710 mV, 2.2252 V and 4.2659 V respectively (about a
1% level difference) with no rejected Newton solutions. The split is made at
Q2's low-impedance emitter-follower output; the real 68 kOhm/6.8 kOhm input arms
remain explicit in the tone solve. At small signal, the actual component network
spans about 9.2x at 80 Hz across the Bass control and 6.2x at 8 kHz across the
Treble control.

The Figure 11-8 voltage and tone stages use trapezoidal capacitor companions.
Backward Euler had added about 2.1 dB of numerical loss at 5--8 kHz at the
96 kHz circuit rate. The complete real-time preamp is now -2.00 dB at 5 kHz and
-8.78 dB at 8 kHz relative to 1 kHz, versus -1.92 and -8.38 dB in ngspice;
doubling and quadrupling the solve rate converge to -1.95/-8.48 and
-1.93/-8.40 dB. The stiff Figure 11-9 transformer/feedback loop deliberately
retains backward Euler: pure trapezoidal companions excited a non-physical
sample-to-sample mode there, which the direct power-circuit tests reject.

Those electrical curves are strongly bunched near the linear pots' end stops.
The module keeps the physical centre and endpoints but maps the Bass panel
position through a signed quarter-power law and Treble through a signed 0.70
power law. This gives every quarter of travel a clearly measurable spectral
change rather than using most of the knob to reach the last few physical
degrees. The historical Treble network supplies only about 4 dB of boost from
its centred position; above centre, the module also parallels its 68 kOhm input
arm progressively down to about 32 kOhm. This is an explicit extended-range
sound-design modification inside the same active-feedback topology, not a
production Peterson component claim or a post-circuit shelf.

The historical 100 kOhm volume loading is present at Q5. The subsequent
post-volume vibrato-feed buffer and lamp/LDR network remain reductions: a
smoothed lamp and equal-power two-channel router drive the independent power
modules. This is the stated boundary of the current Figure 11-8 audio model;
measurements of an original unit are still needed to validate the selected
2N3392 lot and source/load tolerances.

`Drive` is an explicit extension, not a historical Peterson control. It applies
up to 24 dB before the Figure 11-8 input circuit and applies the exact reciprocal
at the schematic's 100 kOhm volume node, before Figure 11-9. This cancels the
knob's linear loudness change while preserving preamp gain loss and harmonics;
it does not use the power amplifier as part of the level compensator. `Output`
is strictly after the nonlinear amplifier. There is no Drive waveshaper or
harmonic generator. An eased-quadratic dB taper keeps the default and moderate
chord range clean, then reaches the deliberate shared-circuit overload range
toward the end stop without concentrating all audible change in its last few
degrees.

The explicit calibration is 0.13 V per model input unit and the module feeds
five times the summed pickups. A hard middle-C peak around 0.0404 model units
therefore presents about 0.42 V peak to the preamp at maximum Drive. A 0.315
model-unit near-pickup/chord stress presents about 0.65 V peak. After the
volume-node reciprocal those cases ask only about 1.5 V and 2.3 V peak,
respectively, from the ideal 57x power stage rather than approaching its 35 V
rail. The complete centred preamp's
6.323 small-signal gain is normalized before Figure 11-9 so existing voltage
calibration is preserved, but nonlinear gain loss is not cancelled. These
voltage-domain mappings are **calibrated**, not measurements of an original
pickup rail.

### Polyphony and distortion placement

The physical pickup rail and Peterson electronics are shared: once note signals
have been summed, any nonlinear preamp or power stage necessarily produces
cross-note intermodulation and sees larger chord peaks. That is different from
several software EP designs. The open-source mda EPiano, for example, applies
its quadratic `overdrive` independently inside the active-voice loop and only
then sums the voices. Its chord remains clearer by construction, but that
placement is not the Suitcase circuit topology. Sampled instruments can likewise
bake bark into each velocity-layer recording before their runtime voice sum.

The DAFx-17 physical model and Falaize/Helie port-Hamiltonian work both locate a
major part of Rhodes single-note spectral enrichment in the nonlinear magnetic
pickup/tine interaction; neither is evidence that a shared Suitcase power stage
should be driven to its rails during ordinary chords. The implementation target
is therefore: per-note fatness and velocity bark primarily in
`ElectricPianoVoice`'s magnetic pickup, realistic headroom in the shared
Peterson path, and explicitly creative chord compression only toward the upper
end of the non-historical `Drive` extension.

The Q1/Q2 and active-tone solves each use an analytic predictor/correction in
the tested range; junction-limited samples continue through checked refinement.
The power-stage solve now always rechecks KCL after its first correction: the
old unconditional acceptance left an absolute error that was small near the
rails but disproportionate on quiet high notes. The regression suite records
zero rejected solutions from small signal through calibrated overload.

## Real-time cost and modulation readiness

`tests/electric_piano_benchmark.cpp` measures one rendered second of a fixed
voice, sixteen fixed voices, sixteen voices plus the shared amplifier, legacy
continuous pitch CV, and active exponential-FM, linear-FM, and phase-modulation
paths. It is a repeatable development microbenchmark rather than a Rack
CPU-meter claim.

The four-times-oversampled pickup remains intact. A process-wide pole-field
table replaces each voice's repeated `tanh` and fractional-power field
evaluations with bilinear interpolation. On the release-Clang
benchmark used during this pass, one fixed voice takes about 6.6 ms, sixteen
fixed voices 106 ms, and the complete sixteen-voice/shared-amplifier path
242 ms per rendered second. Continuously modulated V/oct takes about 233 ms;
active EXP FM, linear TZ FM and PM take approximately 271, 256 and 253 ms.
The 64x physical contact solve contributes only at note onset; held-note cost
is still dominated by pickup and amplifier processing. These figures
are comparative only and will vary with compiler and processor.

The benchmark also contains explicit repeated-strike rows so contact work is
not hidden by a held-note measurement. Reusing the modal contact-point state
between adjacent 64x substeps and caching invariant velocity and inverse-mass
projections leaves render checksums unchanged. On the recorded MinGW/GCC 16.2
run, one voice at 100 strikes/s changed from 86.09 to a median 81.68 ms per
rendered second (-5.1%); sixteen voices at 20 strikes/s each changed from
263.62 to 256.98 ms (-2.5%). Held-voice checksums and cost were unchanged.

The benchmark also isolates the shared amplifier with Vibrato inactive and
active (about 118 and 165 ms respectively on that build). Static Bass, Treble
and Vibrato Speed control laws are cached, values shared by both 2x circuit
frames are evaluated once at host rate, and the lamp/LDR state continues to run
at zero Intensity so enabling Vibrato preserves its physical phase and envelope.
The isolated benchmark precomputes its stimulus outside the timed region. None
of these changes alters a steady-state render checksum. At the Rack boundary,
disconnected modulation and direct
outputs also avoid their otherwise redundant per-voice routing work. The
transistor solves, their iteration criteria, and both oversampling factors are
unchanged.

Adding the separate tone-bar submode initially pushed the three audio-rate
modulation cases to roughly 370/355/352 ms because an eleven-element hot loop
crossed a compiler vectorization boundary. The submode is physical and is not a
creative FM/PM target, so its 4x advance/readout is handled separately while the
existing ten-coordinate modulation loop retains its vectorized layout. This is
numerically identical to the unsplit eleven-mode render.

Pitch is deliberately not cached per strike. A pitch-only update recomputes
modal angles, band limits and the state-preserving velocity transform every
sample, while reusing pitch-invariant coupled-fork eigenvectors and decay radii.
This reduced the benchmark's sixteen-voice continuously V/oct-modulated case
from about 382 ms to about 203 ms without freezing or control-rating pitch.

### Modulation architecture

The Rack boundary uses the following deliberate routing rule:

- V/oct and the single pitch-modulation input use Rack's poly-voltage read. A
  one-channel cable is therefore broadcast to all active notes, while a
  polyphonic cable addresses corresponding voices. Preserved release tails keep
  following a one-channel global source. A stolen polyphonic channel belongs to
  its replacement note, so the old tail instead returns through the existing 4x
  interpolation to zero modulation; feeding it the reassigned voltage would
  retune the released note with the new note's expression. The explicit pedal input
  preserves the model's separate key-held and damper-held states for raw modular
  gates. Rack's MIDI-to-CV module has no pedal output, so the smoke patch adds
  Core MIDI CC-to-CV with cell 1 mapped to unsmoothed CC64. MIDI-to-CV itself
  still holds note gates under CC64 and therefore cannot expose physical key-up
  timing during sustain; this is a limitation of that stock output topology.
- Rack MIDI-to-CV's polyphonic retrigger output feeds the explicit retrigger
  input. A rising retrigger on a still-high channel preserves a differently
  pitched old voice as a tail and strikes the replacement voice. A same-pitch
  retrigger adds a new physical collision to the existing resonator motion.
- The three-way switch routes the modulation input exclusively to exponential
  FM, linear through-zero FM, or phase modulation. The bipolar Depth control
  reverses modulation polarity without changing mode.
- The earlier general-purpose target bay was removed after controlled LFO
  renders showed that most physical parameters were either useful only as
  deliberate static model settings or audibly ineffective under continuous CV.
  Keeping those jacks made the interface larger without creating compelling
  musical behavior.

At full positive depth, exponential FM uses 0.2 octaves/V. Linear FM
uses an additive deviation of 20% of the current positive modal frequency per
volt; consequently -5 V at full positive depth reaches exactly zero on every
key, and more-negative voltage reverses direction. Phase modulation uses 36
degrees/V. Inputs are bounded before they reach the model, and parameter sums
are finite at the Rack boundary.

EXP FM, linear through-zero FM and PM are explicitly **creative electronic
extensions**, not claims that a tine has negative physical resonance. V/oct is
the physical retuning input. The Hunt--Crossley collision, coupled tine/tone-bar
system, damping and positive modal frequencies therefore continue unchanged
under the three creative inputs, while key geometry, graduated hammer
properties and modal mass remain tied to the struck V/oct key.

The three near-fundamental coupled body coordinates receive full modulation.
The seven deliberately inharmonic, short-lived impact coordinates receive none:
audio-rate modulation of those ratios created a dense dissonant sideband cloud,
whereas the nonlinear magnetic pickup already generates a coherent upper
spectrum from modulated body motion. This is a deliberate sound-design routing
choice, not a physical claim about the tine.

No FM topology can make an arbitrary free-running LFO harmonic: sidebands at
`carrier +/- n * modulator` lie on one harmonic grid only when the two
frequencies are commensurate. A controlled 256 Hz render repeated essentially
periodically with a 128 Hz modulator; at 100 Hz, both a plain sine carrier and
the EP were non-periodic. The EP retains a little more roughness because its
coupled body is intentionally split around ratios 0.995, 1.000 and 1.002 and
then passes through the nonlinear pickup. A truly single-carrier synthesizer
mode would be a separate creative model, not a correction to modal FM.

Physical modal coordinates are captured at four evenly spaced pickup-rate
frames (the contact solve remains at 64x). EXP and linear control signals are
reconstructed at 4x and integrated into an independent unit-complex phase
accumulator for each body mode. PM is reconstructed at the same rate and its
frame derivative contributes to pickup velocity, maintaining the Faraday-law
relationship between rendered displacement and EMF. At -100% linear deviation
the body pickup velocity reaches zero continuously; below it the body phase
reverses. The unmodulated impact modes remain audible at that crossing.

Independent accumulators remove the former defect in which one wrapped phase
was multiplied by non-integer modal ratios, making several coordinates jump at
every wrap. Normal increments use a bounded seventh/sixth-order complex rotation
and periodically renormalize; unusually large one-frame changes fall back to
library sine/cosine. Only control types that have actually become nonzero run a
4x control interpolator.

A never-modulated voice retains the established aggregate pickup interpolators,
with no per-mode 4x or trigonometric work. A connected input carrying exactly
0 V is therefore equivalent to no cable. If modulation first becomes nonzero
on an already-ringing note, a 6 ms transition joins the aggregate and per-mode
readouts without a timing step. Once frequency modulation has accumulated
phase, that phase remains when its instantaneous CV returns to zero; abruptly
discarding it would click. A new strike on a fully silent voice resets that
history safely. Any DC PM already present is installed as its initial offset
rather than misread as a one-sample frequency impulse.

## Figure 11-9 power modules and shared supply

`src/models/PetersonPowerAmplifier.hpp` is a coupled ten-unknown nonlinear
circuit solve, not a fitted transfer function. Each stereo channel contains:

- the 6.4 uF input coupling capacitor and Q1 2N3391A error stage;
- the 100 Ohm/5.6 kOhm global feedback network and 1200 pF parallel capacitor;
- Q2, its 47 Ohm emitter path, 22 kOhm shunt, 100 pF collector/base capacitor,
  42 kOhm core-loss path, and the class-A transformer primary;
- both matched 120725/DTG110B **PNP germanium** output devices in the asymmetric
  common-emitter/common-collector topology actually drawn in Figure 11-9;
- both 2.7 Ohm secondary windings, 820 Ohm base/emitter shunts, 0.01 uF
  capacitors, 0.5 Ohm rail paths, the 270 Ohm bleeder, and the electrical
  speaker-load current in the same KCL system.

Ebers-Moll terminal currents with finite forward and reverse beta participate
in an analytic-Jacobian Newton solve. Backward-Euler capacitor companions avoid
parasitic trapezoidal ringing in this stiff feedback circuit. Transistor-
junction-aware step limiting supplies the convergence function normally handled
by a SPICE `pnjlim`; every sample verifies its corrected KCL residual, while only
the steep overload knee needs further linear solves. No `tanh`, `asinh`, polynomial,
harmonic rejection, or post-hoc clip curve creates the amplifier distortion.

The transformer flux is a dynamic state driven by primary voltage. Its
magnetising branch is linear, and secondary current is reflected into Q2
collector KCL. Turns ratio, winding resistance, magnetising inductance, core
loss and reflected loading therefore remain in the circuit; only the former
unmeasured saturation curve has been removed. The service documents contain no
winding/core measurements, so the 1.55 turns ratio and 8 H magnetising
inductance remain explicit **provisional** values. A nonlinear core should only
return after identification from an original T1 demonstrates an audible effect.

### Power-stage reduction audit

The real-time solver has been profiled and compared with controlled ngspice
ablations rather than reduced by deleting components on sight. Earlier profiling
on the development machine took one active mono 48 kHz path from roughly
10.2%/15.2% of one core at default/maximum Drive to 4.8%/9.9%, and direct
0.20/0.52 V peak benchmarks to about 3.9%/6.4% (figures vary with processor and
build). Those figures predate the mandatory post-correction KCL check and are
retained only as optimization history; the complete module needs rebenchmarking.
The savings came from preserving the identified circuit equations while
reducing their evaluation cost:

- reverse-biased junction exponentials below -20 thermal volts are omitted;
  even for the germanium proxy their missed current is below the nonlinear KCL
  tolerance;
- the audio-rate device exponentials use the existing fifth-order exp2
  approximation, whose measured relative error is about 6e-6 over the relevant
  range, while the offline/DC reference retains `std::exp`;
- the unidentified `sinh` magnetisation curve and its `cosh` Jacobian slope are
  absent from the production path; the flux state now stamps a constant 8 H
  magnetising slope;
- device currents and slopes are shared by the residual and analytic Jacobian;
  clean samples normally require one linear correction plus its residual check.

The dense ten-unknown linear solve is now roughly half of the remaining clean
cost. A no-pivot version saved only about 7% of that solve and was rejected
because it weakens numerical robustness for little total benefit. A worthwhile
next solver reduction would eliminate the circuit's linear internal nodes by an
algebraically verified Schur complement or generated sparse factorisation; it
must produce the same terminal equations and overload trajectory.

At 250 Hz the offline ablations give the following sensitivity figures. Values
are output THD for 0.2 V peak / 0.5 V peak at the power-module input:

| Circuit variant | THD at 0.2 V | THD at 0.5 V | Interpretation |
|---|---:|---:|---|
| Full production circuit (linear core) | 0.218% | 1.186% | Reference |
| Former provisional nonlinear core | 0.360% | 6.109% | Unmeasured core curve dominated pushed coloration |
| Resistive speaker load | 0.175% | 0.914% | Load changes overload balance, not clean level dramatically |
| No 0.01 uF output-device capacitors | 0.218% | 1.186% | Negligible at 250 Hz, retained for HF behaviour/stability |

The former provisional magnetisation curve raised 0.5 V THD by more than five
times and introduced a large DC/asymmetric shift without measurement support.
Production therefore retains the electrically significant transformer states
but keeps the core linear. The reactive load and output compensation capacitors
cost little enough to retain.
The ablation mode and its raw harmonic output are reproducible with
`uv run --with numpy python tests/python/reference_peterson_power_spice.py
--ablations` when ngspice is installed.

The real-time implementation is checked against
`tests/python/reference_peterson_power_spice.py`, a durable ngspice netlist of
the same topology. Its DTG110B Gummel-Poon proxy is constrained by the scanned
short-form data (PNP germanium, hFE 65--300 at 1 A, fT about 500 kHz), but no
manufacturer SPICE model has been found. The real-time Ebers-Moll saturation
current and effective thermal slope are reduced parameters fitted to the
low-current and power-region curve of that fuller proxy, rather than copied as
if they were manufacturer values. At 250 Hz and 0.05 V peak, ngspice gives
1.801 V AC RMS with H2/H3/H4 of 0.138%/0.096%/0.076%; real time gives 1.919 V
with 0.148%/0.096%/0.072%. At 0.5 V peak, ngspice gives 19.129 V AC RMS with
H2/H3 of 0.86%/0.81%; real time gives 18.418 V with 4.28%/3.58%. Clean
crossover behaviour is therefore tightly constrained, while pushed-device
harmonics retain wider tolerance because no measured 120725 model exists.

Both modules draw from a common +/-35 V supply with finite droop and recovery.
Figure 11-10's 3000 uF reservoirs support the topology, while effective winding
resistance and recharge timing remain **provisional**. Rack cable protection is
applied only after decimation and Drive compensation; it cannot cause or shape
normal amplifier overload.

## Electrical reactive speaker load

The service schematic specifies a 16 Ohm speaker load per power module but no
impedance curve. A resistor alone cannot produce the frequency-dependent current
and back-EMF seen by the output pair. Each channel therefore uses the standard
linear moving-coil equations

`v = Re i + Le di/dt + Bl u`

`Mms du/dt = Bl i - Rms u - x/Cms`, with `dx/dt = u`.

The three states are trapezoidally discretized at the 2x amplifier rate. The corresponding
small-signal terminal impedance is

`Z(w) = Re + jwLe + Bl^2 / (Rms + jwMms + 1/(jwCms))`.

Current provisional values are `Re=12.8 Ohm`, `Le=0.55 mH`, `Mms=18 g`,
`Cms=0.2502 mm/N`, `Rms=2.12 N s/m`, and `Bl=10.5 T m`. They produce a
16-Ohm-class DC region, a roughly 75 Hz motional peak around 65 Ohm, a return
toward nominal impedance in the midband, and an inductive high-frequency rise.

Suhr and Tone King reactive-load documentation is relevant here only as evidence
that a useful amplifier load must reproduce the speaker's resonant and
inductive impedance, separately from any IR/cabinet filter. Their products are
8-Ohm guitar-cabinet targets driven mainly by tube amplifiers; their exact curve
and the strength of their amplifier/load interaction are **not** Rhodes data.
The globally fed-back solid-state Peterson module should show much less clean
voltage-response coloration from the load than a tube output transformer. Here
the load matters chiefly through output-device current, junction drop, rail
headroom and supply recovery; it must not become a disguised cabinet EQ.
The current Rhodes parameters are therefore explicitly **provisional** until an
original Suitcase speaker pair is measured in its enclosure. Useful measurements
are magnitude and phase from 10 Hz to 30 kHz at small signal, plus compression
and inductance modulation at several drive levels.

## Validation currently in the DSP suite

Tests cover contact stability and natural separation across sample rates;
finite contact-footprint bounds, velocity independence and ultrasonic
contact-only modal loading; coupled-fork mass orthogonality,
hybridization and support reaction; published F1/F3 ratio and decay anchors;
spring-loaded eigenfunction contact residues; separate reciprocal tone-bar
submode frequency, decay and displacement anchors at F1/F3;
the short-tine shear correction and oversampled modal-band taper; the symmetric
73-key Strike law; keyboard level/decay span; four-key lower-register level,
envelope, H2/H3 and velocity-bark fingerprints; the 73-key
multi-velocity pickup fit, direct-field/LUT agreement and register-wide
Proximity level window; a C6/C7/top-key harmonic-bed guard against isolated
inharmonic lines; monotonic velocity energy and spectral growth; pickup harmonic
and intermodulation growth; audible Tone, Hammer, Bell, midpoint-to-endpoint
Coupling and Decay ranges;
polyphonic mechanical-event timing; explicit circuit/Rack voltage conversion;
the full Drive/THD trajectory; an actual eight-note maximum-Drive render that
guards power-rail headroom and chord crest; direct Figure 11-8 input and
active-tone preamp clean/overload level, THD and convergence against ngspice;
Figure 11-8 5/8 kHz response agreement with ngspice;
amplifier clean/default/
driven distortion; reactive-load DC, resonance, midband and inductive regions;
direct real-time/ngspice level and harmonic agreement; zero nonlinear
solver failures through rail overload; stereo optical routing; shared-supply
recovery; bounded rails; 7 kHz maximum-Drive nonlinear-alias rejection; a
default-Drive 4.2 kHz/C8-region folded-crossover check; and 48/96 kHz
cross-rate residual checks on a rendered single note and four-note chord at
maximum Drive; independent FM phase-wrap and slow-PM continuity; a true
linear-through-zero body crossing; and a 48/96 kHz audio-rate harmonic-FM
residual comparison.

These are invariants and regression checks, not proof of similarity to a
particular piano. The next report should add source recordings, measurement
conditions, plots, tolerances, parameter fits, and a versioned render corpus.

## References

1. Fender Rhodes, [Service Manual](https://www.fenderrhodes.com/service/manual.html), especially Chapters 1, 4, 7 and Figures 11-8 to 11-10.
2. K. H. Hunt and F. R. E. Crossley, [“Coefficient of Restitution Interpreted as Damping in Vibroimpact”](https://doi.org/10.1115/1.3423596), 1975.
3. A. Falaize and T. Hélie, [“Passive simulation of the nonlinear port-Hamiltonian modeling of a Rhodes Piano”](https://doi.org/10.1016/j.jsv.2016.11.008), 2017.
4. F. Pfeifle, [“Real-time Physical Model of a Wurlitzer and Rhodes Electric Piano”](https://www.dafx17.eca.ed.ac.uk/papers/DAFx17_paper_79.pdf), DAFx-17.
5. Vintage Vibe, [Peterson four-pin preamplifier rebuild documentation](https://www.vintagevibe.com/products/fender-rhodes-4-pin-pre-amp-rebuild-kit), a higher-resolution annotated copy of Figure 11-8 used to disambiguate component values and transistor pin orientation.
6. L. Gabrielli et al., [“The Rhodes electric piano: Analysis and simulation of the inharmonic overtones”](https://doi.org/10.1121/10.0002002), 2020.
7. D. T. Yeh and J. O. Smith, [“Simulating Guitar Distortion Circuits Using Wave Digital and Nonlinear State-Space Formulations”](https://www.dafx.de/paper-archive/2008/papers/dafx08_04.pdf), DAFx-08.
8. Paul Kellett, [mda EPiano open-source implementation](https://sources.debian.org/src/mda-lv2/1.2.6-1/src/mdaEPiano.cpp), especially the per-voice overdrive at lines 971-989, used as a polyphonic architecture comparison rather than a circuit reference.
9. Rhodes Music, [Vari-Amp](https://rhodesmusic.com/rhodes-vari-amp-plugin/) and its [manual](https://rhodesmusic.com/wp-content/uploads/2025/08/Vari-Amp-User-Manual.pdf), showing separate vintage-preamp, drive-type, parallel-mix, and amp sections in a modern commercial design.
10. Suhr, [Reactive Load product description](https://www.suhr.com/electronics/tone-tools/suhr-reactive-load/) and [manual](https://www.suhr.com/wp-content/uploads/2020/07/Reactive-Load-070620.pdf).
11. Tone King, [Ironman II manual](https://www.toneking.com/wp-content/uploads/2019/09/Ironman-II-Manual-.pdf), reactive-load discussion and example impedance curves.
12. Celestion, [Guitar Loudspeaker Catalogue](https://celestion.com/wp-content/uploads/2020/03/Guitar_Speaker_Catalogue.pdf), nominal impedances and representative driver data.
13. OpenWurli, [circuit-first open implementation](https://github.com/hal0zer0/openwurli), used as an independent architecture comparison, not as Rhodes measurement data.
14. General Transistor, *DTG110B Germanium PNP High Power Transistor* short-form data sheet (scan retained during model development); hFE and fT constrain the provisional device model.
15. T. L. Sonderbo, [*Real-time physical model of the Rhodes Electric Piano*](https://projekter.aau.dk/projekter/files/719175589/Master_Thesis_Real_time_physical_model_of_the_Rhodes_Electric_Piano_Tobias_Sonderbo.pdf), Aalborg University, 2024; tables 3.1--3.3 and sections 3.2--3.4 supply the tine, Hunt--Crossley, and damped-spring quantities used in SI units.
16. Steve Woodyard, [Rhodes striking-line service notes](https://www.fenderrhodes.com/service/escapement.html) and [tine/pickup settings](https://www.fenderrhodes.com/service/tine-settings.html), used to distinguish a clear maximum-power strike from taper/tip “thunk” and to cross-check practical pickup gaps.
17. A. Falaize, [Rhodes model article and downloadable companion code](https://afalaize.github.io/posts/rhodes/), used to verify the published A4 contact geometry and the 15 mm hammer-width integration implemented by the reference model.
18. John Learman, [jRhodes3d direct-harp sample corpus](https://github.com/jlearman/jRhodes3d-wav), 1977 Mark I Stage 73 recordings used for aligned register/velocity spectrogram trends. The files are mildly EQ'd and noise-reduced, so absolute spectral tilt is not treated as raw SLDV data.
