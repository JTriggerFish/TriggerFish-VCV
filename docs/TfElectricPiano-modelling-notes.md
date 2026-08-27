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
| Hammer contact while engaged | 16x host | nonlinear physical model |
| Magnetic pickup | 4x host | reduced field geometry |
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

`Vel Curve` applies a controller-response curve. `Dynamics` maps the result to
the physical incoming hammer-speed range. A very small positive velocity is
treated as controller residue before a collision starts; an active collision
is never terminated by a timer. **Calibrated:** the default Dynamics value is
the full input range. The present action is not a multi-lever key escapement
simulation, so hammer speed as a function of key travel remains an open model.

## Hammer contact

The finite hammer mass and all participating modal velocities are coupled in
both directions. For compression `q > 0`, the contact force has the
Hunt-Crossley form

`F = k q^n (1 + alpha qDot)`, clamped to a non-tensile force.

The hardness control changes `k`, `n`, loss, and hammer mass rather than adding
a brightness EQ. Separation follows from zero compression/receding contact.
**Source:** Hunt and Crossley's nonlinear damping term and the Rhodes hammer-
tine topology. **Calibrated:** neoprene stiffness, exponent, loss and keyboard
scaling. These require force/displacement or contact-duration measurements from
tips of known hardness for absolute validation.

## Tine, tone bar and mounting block

The fundamental is a two-coordinate asymmetric fork. Tine and tone-bar mass,
stiffness and common-base coupling form a 2x2 generalized eigenproblem. Its
normal frequencies, modal masses, strike projections, tine pickup projections,
and support-reaction loss follow from the same eigenvectors. Coupling therefore
changes energy exchange and decay rather than crossfading two oscillators.

Higher signed attack modes use cantilever-like ratios, strike-position mode
shapes, and keyboard-dependent damping. The inharmonic ladder is retained only
while below the anti-alias taper. **Source:** the Rhodes service manual's
unequal-leg tuning fork and published modal studies. **Derived:** the coupled
eigenproblem and mass orthogonality. **Calibrated/provisional:** absolute modal
masses, base stiffness, higher-mode amplitudes, and loss curves. Laser-vibrometer
data across several keys would supersede these calibrations.

## Magnetic pickup

Two transverse tine coordinates move through an asymmetric, two-protuberance
magnetic-pole approximation. Flux is a nonlinear function of the instantaneous
tine position and the generated voltage is its time derivative, consistent
with Faraday's law. `Tone` moves the vertical resting alignment; `Proximity`
changes the pole gap. A partial small-signal normalization stops Proximity from
becoming merely another level control while preserving curvature and overload.
Across its full range the current calibration changes isolated-note energy by
about 1.64x (roughly 2.2 dB in amplitude), rather than the earlier 11--13 dB.

**Source:** published Rhodes field/geometry models and the service adjustment
procedure. **Derived:** flux differentiation and velocity-dependent harmonic
and intermodulation generation. **Calibrated:** pole dimensions, scale and
per-key sensitivity. A measured two-dimensional flux map and pickup impedance
are still needed. The model does not contain a separate bark waveshaper.

## Damper and mechanical noise

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
  output capacitor are backward-Euler companions in the same Newton solve;
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
voice, sixteen fixed voices, sixteen voices plus the shared amplifier, and
sixteen continuously pitch-modulated voices. It is a repeatable development
microbenchmark rather than a Rack CPU-meter claim.

The four-times-oversampled pickup remains intact. Its two magnetic edges share
one radial falloff calculation, and the remaining `tanh` and fractional-power
evaluations use bounded approximations with explicit regression limits of
`1.1e-4` absolute and `5.7e-5` relative respectively. On the release-Clang
benchmark used during this pass, sixteen fixed voices fell from about 213 ms to
110 ms per rendered second. Adding the shared nonlinear amplifier currently
takes the complete sixteen-voice path to about 247 ms on that same machine.
These figures are comparative only and will vary with compiler and processor.

Pitch is deliberately not cached per strike. A pitch-only update recomputes
modal angles, band limits and the state-preserving velocity transform every
sample, while reusing pitch-invariant coupled-fork eigenvectors and decay radii.
This reduced the benchmark's sixteen-voice continuously modulated case from
about 382 ms to 189 ms without freezing or control-rating pitch. Future pitch
CV/FM can therefore remain sample-accurate; other modulation targets should be
classified similarly by whether they affect oscillator angle, decay, coupling,
or pickup geometry.

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
coupled-fork mass orthogonality, hybridization and support reaction; keyboard
level/decay span; monotonic velocity energy and spectral growth; pickup harmonic
and intermodulation growth; audible Tone, Hammer, Coupling and Decay ranges;
polyphonic mechanical-event timing; explicit circuit/Rack voltage conversion;
the full Drive/THD trajectory; an actual eight-note maximum-Drive render that
guards power-rail headroom and chord crest; direct Figure 11-8 input and
active-tone preamp clean/overload level, THD and convergence against ngspice;
amplifier clean/default/
driven distortion; reactive-load DC, resonance, midband and inductive regions;
direct real-time/ngspice level and harmonic agreement; zero nonlinear
solver failures through rail overload; stereo optical routing; shared-supply
recovery; bounded rails; 7 kHz maximum-Drive nonlinear-alias rejection; a
default-Drive 4.2 kHz/C8-region folded-crossover check; and 48/96 kHz
cross-rate residual checks on a rendered single note and four-note chord at
maximum Drive.

These are invariants and regression checks, not proof of similarity to a
particular piano. The next report should add source recordings, measurement
conditions, plots, tolerances, parameter fits, and a versioned render corpus.

## References

1. Fender Rhodes, [Service Manual](https://www.fenderrhodes.com/service/manual.html), especially Chapters 1, 4, 7 and Figures 11-8 to 11-10.
2. K. H. Hunt and F. R. E. Crossley, [“Coefficient of Restitution Interpreted as Damping in Vibroimpact”](https://doi.org/10.1115/1.3423596), 1975.
3. A. Falaize and T. Hélie, [“Passive simulation of the nonlinear port-Hamiltonian modeling of a Rhodes Piano”](https://doi.org/10.1016/j.jsv.2016.11.008), 2017.
4. J. Rauhala et al., [“Real-time Physical Model of a Wurlitzer and Rhodes Electric Piano”](https://www.dafx17.eca.ed.ac.uk/papers/DAFx17_paper_79.pdf), DAFx-17.
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
