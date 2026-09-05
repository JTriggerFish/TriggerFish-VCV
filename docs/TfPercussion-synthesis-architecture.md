# Modular percussion synthesis and calibration architecture

## Purpose

TriggerFish should grow into a modular perceptual percussion synthesizer rather
than a collection of unrelated drum voices. Instruments should be assembled
from deterministic, independently tested DSP components, while shared analysis
and calibration tools fit those assemblies to real recordings.

The target is not an exact finite-element simulation. It is a controllable
model whose contact, energy transfer, resonance, radiation, and loss correspond
to audible parts of a real instrument. A ride, hi-hat, snare, kick, tom, or
other struck object can then reuse the same components with different topology
and fitted parameters.

The previous cymbal engines and their fitted parameters have been removed.
The Rack modules remain silent interface shells while the new components are
implemented and validated. No deleted model is an A/B or calibration baseline;
comparisons are always against real reference recordings.

Current clean-slate component coverage and its remaining acceptance gates are
tracked in [the DSP implementation status](TfPercussion-dsp-implementation-status.md).

Development starts with a crash. Its high-energy swell is the harder single-
cymbal case and the resulting body can be reused for strong ride articulations.
The first same-object corpus has five strike articulations, 9--16 velocity
layers, and multiple round robins.

## Code organization contract

The replacement must remain readable as it grows:

- reusable DSP lives under `src/tfdsp/percussion/`, grouped by primitive rather
  than by instrument;
- one focused public component per header/source pair, with implementation
  details kept private;
- ride, hi-hat, snare, and kick files only compose components and map controls;
  they do not contain filter, oscillator, delay, or fitting algorithms;
- tests mirror components, while graph and Rack integration tests remain
  separate;
- analysis separates audio I/O, datasets, segmentation, descriptors,
  objectives, optimization, reporting, and audition assembly;
- shared helpers must express a real invariant or unit, not merely shorten one
  caller; and
- names include units where ambiguity is possible (`frequencyHz`,
  `durationSeconds`, `delaySamples`, `gainLinear`).

As a review trigger, a component file approaching roughly 300 lines or a
function approaching roughly 50 lines should be split unless keeping it whole
makes the invariant materially clearer. This is not a mechanical limit, but a
large exception must be justified in the file. Public APIs should stay smaller
than their implementations and avoid exposing mutable internal state.

The word *graph* in this document describes signal flow and the serialized fit
topology. It does not require a runtime graph engine. Production instruments use
statically typed component members and direct, inlinable C++ function calls in
a fixed order. The audio thread performs no node traversal, virtual dispatch,
allocation, or topology discovery. Test/offline variants may compose the same
components differently, while production topology choices are compiled or
prepared outside audio processing.

## Common event and control model

Every instrument accepts a common contact event containing at least:

- strength or velocity;
- contact location;
- hardness;
- implement type or a continuous implement coordinate;
- finite contact duration and, where relevant, the duration of a held gate.

Persistent controls such as location, hardness, implement, mute, pedal, size,
and tuning have panel defaults that CV can override. Hit and strength remain
event inputs. A trigger and a gate are both valid: a trigger creates one finite
contact, while a gate may also describe a continuing gesture such as a brush
sweep. Mute is passive and must never inject energy.

## Reusable contact excitation

One contact model combines three primitives. Their mix is selected by the
implement model rather than added unconditionally.

1. **Finite force pulse.** A normalized half-sine
   compression/release pulse represents coherent contact. It is particularly
   useful for sticks and compliant mallets.
2. **Tonal contact chirp.** A very short damped oscillator with a pitch or FM
   envelope supplies the hard-tip tick and focused bell articulation.
3. **Distributed micro-contact process.** Overlapping noise or grains represent
   many small contacts. They are contact statistics, not audible granular
   events. A compact cluster can describe a brush tap; a continuing stream can
   describe a brush sweep.

Typical implement mappings are:

| Implement | Finite pulse | Tonal chirp | Micro-contacts |
| --- | --- | --- | --- |
| Hard stick | short and strong | strong | minimal |
| Soft stick | broader | subtle | low |
| Mallet | broad half-sine | little or none | dense, smooth, short |
| Brush tap | weak | none | compact cluster |
| Brush sweep | minimal | none | extended correlated stream |

Hardness changes contact duration, chirp strength, and contact bandwidth.
Velocity changes incident energy and may move a nonlinear body into a different
operating regime. Location changes coupling into the body; it does not turn one
implement into another.

The contact component produces two explicit routes:

```text
contact event -> direct contact radiation
              -> body drive
```

This permits an audible stick transient while ensuring that noise used to
excite a body is not mistaken for an independently mixed noise layer.

## Detailed architecture documents

The remaining design is split by responsibility:

- [Nonlinear resonator architecture](TfPercussion-nonlinear-resonator-architecture.md):
  self-contained signal flow, contact and bloom behaviour, the pooled 512-state
  unified stochastic modal field, energy invariants, controls, and current
  limitations.
- [Metallic percussion DSP components](TfPercussion-metal-components.md):
  dispersion, feedback causality, resonators, shifters, fractional delay,
  radiation, and live loss.
- [Percussion instrument topologies](TfPercussion-instrument-topologies.md):
  reusable membrane/contact structures and candidate cymbal, hi-hat, snare,
  and kick assemblies.
- [Percussion DSP quality and calibration](TfPercussion-quality-and-calibration.md):
  component acceptance, numerical references, integration tests, rendering,
  and calibration.
- [Crash fitting methodology](TfCrash-fitting-methodology.md): the first-object
  manual fitting baseline and its reproducibility requirements.
- [Percussion ear-fitting workbench](TfPercussion-ear-fitting-workbench.md):
  manual fitting controls, browser architecture, snapshots, and safety.
- [Modular percussion patches](TfPercussion-patch-architecture.md): versioned
  JSON graphs, module ports, workbench lowering, and reduced Rack controls.
- [Percussion analysis toolkit](TfPercussion-analysis-toolkit.md): canonical
  analysis contracts, dataset policy, validation diagnostics, and the
  Plotly/WebAssembly reporting boundary.
- [Cymbal reference corpus](TfRide-reference-corpus.md): source-data inventory and
  provenance.

## Open architectural questions

- Whether the fixed 17-state packet provides enough density at every anchor or
  should become an adaptive allocation after perceptual validation.
- Whether local orthogonal modal exchange is sufficient, and how much exchange
  improves metallic evolution before it blurs useful ridges or causes
  implausible redistribution between decay regions.
- Whether slow delay modulation remains necessary once self-phase distortion
  and stochastic modal diffusion are present.
- Whether the slow and self-modulated reads should ultimately be combined.
- Whether bell/plate and later shell/membrane subsystems need separate
  dispersion loops or only separate resonator projections.
- Whether membrane bodies should use the same sparse-modal/dense-residual
  decomposition selected for cymbals or a compact FM body at each complexity
  level.
- How much of a kick's fitted far-field residual belongs to the instrument
  radiation model and how much belongs to the external room renderer.
- Whether shell and port resonances warrant dedicated components after the
  first coupled two-head fits.

## Cymbal representation decision

The current proposed crash body is one stochastic modal field. Zero to 32
editable anchors reserve coherent centres, while deterministic satellite pairs
share the remaining capacity of one 512-state pool. Direct contact and
intrinsic bloom excite the same stored field.
Normalized centre-to-satellite drive, ERB-domain spread, magnitude-preserving
phase diffusion, and local orthogonal energy exchange provide a continuous
resolved-ridge-to-wash trajectory. Raw dispersion remains excitation and
analysis only. The former sparse modal bank, statistical cloud, and turbulent
reservoir remain a legacy workbench comparator rather than three parts of the
proposed output. The complete graph and equations are in the
[nonlinear resonator architecture](TfPercussion-nonlinear-resonator-architecture.md).

All percussion bodies are mono unless a physical topology explicitly contains
multiple bodies. Mono output remains a complete supported presentation.
Optional stereo is derived at the observation stage from shared mono source
taps, so width and microphone perspective cannot alter synthesis state or
calibration of the object itself.
