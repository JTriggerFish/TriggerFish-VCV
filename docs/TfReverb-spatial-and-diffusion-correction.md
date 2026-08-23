# TfReverb stereo spatial and diffusion correction

## Status

This document specifies the production changes required to make TfReverb a
coherent positioned-source, stereo-output room reverb. It replaces the earlier
proposal. The changes are not yet implemented.

The important product decision is:

> TfReverb renders conventional two-speaker stereo from one coincident
> listener reference point. It is not a binaural renderer or a virtual stereo
> microphone pair.

That decision determines how direct sound, image-source reflections, and the
late field are decoded. It also avoids mixing incompatible amplitude-panning,
spaced-microphone, and Ambisonic models in one signal path.

## Goals and non-goals

The corrected renderer must provide:

- independent positioning of every mono/polyphonic source;
- physically calculated image-source reflection timing and lateral direction;
- a smooth transition from discrete reflections into a spacious late field;
- a late tail that initially reflects source and listener geometry, then loses
  individual-source direction as it becomes diffuse;
- stable, useful interchannel decorrelation without an arbitrary stereo pattern;
- a Diffusion control that changes scattering strength as well as delay span;
- deterministic output for a static scene and click-free parameter automation;
- good mono compatibility and bounded gain throughout the room.

The following remain out of scope:

- HRTFs, head shadow, pinna filtering, and head tracking;
- literal front/rear or elevation localization over two speakers;
- spaced A/B, ORTF, dummy-head, or other virtual microphone rendering;
- surround or Ambisonic output;
- assigning physical directions to FDN state variables that have no directional
  meaning in the current topology.

Front/rear and elevation still affect propagation distance, wall encounters,
material loss, and response timing. They simply do not receive invented stereo
cues that two speakers cannot reproduce reliably.

## Problems in the current implementation

The existing code contains useful geometric information, but the complete
signal path is inconsistent:

1. Each ER path already has a physical direction and equal-power stereo gains,
   but the FIR is placed at absolute propagation time while dry sound is
   immediate.
2. Every dry source is summed into centred mono before the dry/wet mix.
3. The late input and output connections preserve relative wall timing, but
   discard common timing and use normalized directional vectors as though they
   were also a complete distance and coupling model.
4. The late network is decoded through a hand-written `SidePattern` rather than
   through the actual lateral direction of each wall port.
5. `PositionBalanceGains` changes ER and tail in opposite directions while the
   dominant direct signal is unchanged.
6. Diffusion changes only the velvet-delay span; every Hadamard transform remains
   fully mixed at every setting.
7. Existing spatial tests mostly establish that two responses differ. They do
   not establish direction, symmetry, timing, coherence, or convergence toward
   a diffuse field.

The ER stereo pan itself is not missing. The current calculation in
`EnumerateEarlyReflectionPaths()` is the correct starting point for the chosen
speaker-stereo renderer. It needs a precise horizontal convention and stronger
tests, not replacement with two virtual ears.

## Coordinate convention

All source and listener controls remain normalized XYZ positions in the room.
After conversion to metres:

- +X is right;
- -Y is forward, toward the top of the room plan;
- +Z is up;
- the listener has no rotating head orientation;
- stereo lateral direction is calculated from the horizontal XY projection.

For a horizontal vector `v = (x, y)`, define

```text
lateral(v) = x / sqrt(x^2 + y^2)       when x^2 + y^2 > epsilon
             0                         otherwise
```

`lateral` is -1 at hard left, 0 at front or rear, and +1 at hard right.
Elevation must not pull an otherwise lateral source toward the centre; Z is
therefore excluded from this calculation.

The constant-power stereo decoder is

```text
angle = (pi / 4) * (lateral + 1)
gainL = cos(angle)
gainR = sin(angle)
```

so `gainL^2 + gainR^2 = 1` for every decoded source or path.

## Corrected signal path

```text
independent mono sources and XYZ positions
        |
        +-- positioned direct decode ----------------------------> stereo direct
        |
        +-- user wet pre-delay
             |
             +-- zero-latency head/tail image-source FIR
             |      -- direct-relative reflection taps ----------> stereo ER
             |
             +-- source handoff-alignment delay
                    -- six normalized source-wall sends
                    -- 16-line room-scaled FDN and VFM
                    -- six wall-domain outputs
                    -- six physical listener-wall connections
                    -- wall-direction stereo decode -------------> stereo tail

stereo ER + stereo tail -- EARLY/TAIL trims -- wet filters -- WIDTH --> wet
stereo direct + wet -- constant-power MIX -- output LEVEL -----------> output
```

The module remains zero-latency on its direct path. The ER convolver must
therefore use a zero-latency head plus partitioned tail, described below. This
avoids delaying every dry signal by the current 128-sample FFT partition and
avoids non-causal compensation when an ER has less than 128 samples of excess
delay.

## Direct path

### Per-source rendering

The DSP layer must retain every source independently until direct stereo decode.
For source `s`:

1. Convert source and listener positions to metres.
2. Calculate the source-minus-listener horizontal vector.
3. Decode it with the constant-power law above.
4. Apply the bounded direct-distance gain below.
5. Add the stereo pair to the direct accumulator.

The dry API must consequently become stereo. A suitable interface is:

```cpp
struct RoomReverbFrame {
  StereoFrame direct;
  StereoFrame wet;
};

StereoFrame MixReverbOutput(const StereoFrame& direct,
                            const StereoFrame& wet,
                            const ReverbOutputGains& gains);
```

`TfReverb.cpp` must no longer create a scalar `dry` sum.

### Direct mode and patch compatibility

Add a serialized context setting:

```cpp
enum class DirectMode { CentredLegacy, Positioned };
```

- New module instances default to `Positioned`.
- Existing patches that do not contain the setting load as `CentredLegacy`.
- `CentredLegacy` reproduces the current immediate, unattenuated mono sum in
  both channels.
- `Positioned` uses independent panning and bounded distance gain.

The menu and manual must call these modes **Centred dry (legacy)** and
**Positioned direct sound**. A centred mono sum must not be described as
physically positioned output.

### Bounded distance law

Changing Size should not turn the room pad into a large dry-level control.
Distance is therefore expressed relative to the room's characteristic length:

```text
roomScale = cbrt(roomVolume)
rho       = sourceListenerDistance / roomScale
rhoRef    = factoryDistance / cbrt(factoryVolume)

gDirect = clamp(rhoRef / max(rho, 0.05), 0.25, 2.0)
```

The initial limits are -12 dB and +6 dB. The factory source/listener placement
has unity gain. Uniformly scaling the room while preserving normalized source
and listener positions leaves this gain approximately unchanged.

These constants may be adjusted once during validation, but the final values
must be fixed in the implementation and tests. They must not be optimized as
part of the FDN coefficient search.

No interchannel time delay is added to the direct sound. This is conventional
speaker amplitude panning and remains mono-compatible.

## Early reflections

### One listener reference, not two ears

The image-source geometry continues to use one listener position. Each image
source produces one propagation path whose horizontal arrival direction is
decoded with the constant-power law.

Using two receiver points would instead simulate a spaced microphone pair. It
would introduce frequency-dependent interchannel phase, require an explicit
spacing and orientation, and create possible mono-sum comb filtering. Without
head shadow or HRTFs it would still not be binaural. That is a valid future
effect mode, but it is not part of this correction.

### Direct-relative timing

For every path:

```text
dDirect = distance(source, listener)
dPath   = distance(imageSource, listener)
tExcess = max(0, (dPath - dDirect) / speedOfSound)
```

The FIR tap must use `tExcess`, not `propagationSeconds`:

```text
desiredTapSamples = tExcess * sampleRate
```

The current convolver adds one 128-sample partition of latency and compensates
for it by moving sufficiently late FIR taps earlier. That cannot represent an
ER whose excess delay is less than one partition without either delaying the
direct path or making the FIR non-causal.

Replace it with a zero-latency head/tail convolver:

1. The first 128 FIR coefficients for every source/output kernel form the head.
2. Process the head directly from a 128-sample per-source input history.
3. Coefficients from sample 128 onward use the existing uniform partitioned
   convolver, shifted down by 128 coefficient samples so its inherent
   128-sample latency restores their exact requested output time.
4. Sum head and tail results before FIR-bank transition gain is applied.
5. Prepared head and tail banks change through the same numerical crossfade;
   neither input history may be cleared during a scene change.

The worst-case direct head is 128 taps by two outputs by eight active sources.
Benchmark that path at the production polyphony. If it is too expensive, use a
sparse head or non-uniform small first partitions, but retain exact zero-latency
semantics.

`appliedLatencyCompensationSamples` and `ResidualLatencySeconds()` should be
removed from the ER response if they are no longer used elsewhere. Tests and
handoff analysis must operate entirely in direct-relative time.

The user Pre-delay is additional creative wet delay:

```text
audible ER time relative to direct = userPreDelay + tExcess
```

### Direction and gain

For every retained image path:

1. Calculate `pathVector = imagePosition - listenerPosition`.
2. Calculate lateral direction from its XY projection.
3. Calculate constant-power L/R gains.
4. Multiply both channels by the existing distance spreading, air loss, and
   four-band material response.

Do not blend later image paths toward random or synthetic stereo directions.
The transition becomes diffuse through the increasing number and overlap of
physically directed paths. Changing individual directions would weaken room
symmetry and introduce another arbitrary spatial pattern.

### ER handoff taper

The detected per-source handoff remains an interval:

```text
tStart = sourceHandoff.finalStartSeconds
tEnd   = sourceHandoff.finalEndSeconds
```

For a kernel sample at direct-relative time `t`, apply

```text
ER envelope = 1                                      when t <= tStart
              cos(pi/2 * smoothstep((t-tStart) /
                                     (tEnd-tStart))) when tStart < t < tEnd
              0                                      when t >= tEnd
```

The envelope is identical in both channels and therefore cannot move a
reflection while fading it.

## Late-field input and ER/late alignment

### Six wall ports remain the physical interface

The six late ports represent the room's -X, +X, -Y, +Y, -Z, and +Z walls. They
are meaningful geometric signals. The 16 FDN lines are internal mixed state and
must not be labelled as directions.

For each source, retain the existing raw solid-angle proxy and calculate

```text
unitWallWeights = rawWallWeights / norm(rawWallWeights)
relativeDelay[w] = (distanceToWall[w] - minimumDistanceToWall) / c
```

The normalized vector is desirable: it gives a bounded passive injection
direction. Overall source energy is supplied separately by the ER handoff match
`currentTailSend[source]`. Do not apply `PositionBalanceGains` on top of it.

This separation is explicit:

```text
late source send = delayedSource
                 * handoffMatchedScalar
                 * unitWallWeight[w]
```

Normalization is therefore not treated as a distance model, and a second
equal-and-opposite ER/tail distance proxy is unnecessary.

### Source-specific alignment delay

For an impulse injected into the current late topology, define
`tIntrinsicFirst` as the first non-zero wall-domain output time excluding user
Pre-delay. Initially calculate it from the minimum main delay and minimum
connection delays; verify it against a rendered late-only impulse response.

Each source receives a transition delay before its six wall sends:

```text
tAlign[source] = max(0,
                     handoff.finalStartSeconds - tIntrinsicFirst)
```

The audible late onset relative to the logical direct reference is then

```text
userPreDelay + tAlign[source] + tIntrinsicFirst
```

The existing old/new fixed-read-head crossfade mechanism must be used when
`tAlign` changes. It must never clear source, wall, or FDN delay buffers.

If `tIntrinsicFirst` is later than the handoff start for any supported Size,
the test must fail. The correction is then to shorten or restructure the late
network at that Size; a negative delay or abrupt ER extension is not allowed.

The late response must naturally build through `[tStart, tEnd]` while the ER
envelope falls. A global output gate is not suitable because it would also gate
reverberation already stored from earlier notes.

## Late-field output and stereo decode

### Keep the wall-domain output

Continue to project the 16 delayed FDN values back into the same six-dimensional
wall basis:

```text
outputWalls = transpose(WallProjection) * fdnDelayed
```

`WallProjection` has orthonormal columns after its `1/sqrt(16)` scale. Injection
and extraction therefore remain passive. The complete VFM may populate the ten
complementary FDN coordinates, but those coordinates do not acquire physical
directions merely because they exist.

No 12- or 16-source intermediate renderer is added for stereo output. With no
direction-dependent filtering or multichannel reproduction it would collapse
algebraically to another 2-by-16 output matrix and would not increase spatial
resolution.

### Physical six-wall speaker decode

Keep the listener's normalized wall weights and relative wall delays. Replace
`SidePattern` and its Gram-Schmidt construction with constant-power decoding of
the actual wall directions.

For the current wall order, use

```text
wallLateral = {-1, +1, 0, 0, 0, 0}
```

The front, rear, floor, and ceiling ports decode centrally because this product
does not claim front/rear or height localization. For wall `w`:

```text
panL[w], panR[w] = constantPowerPan(wallLateral[w])

decoderL[w] = listenerWallWeight[w] * panL[w]
decoderR[w] = listenerWallWeight[w] * panR[w]
```

Because the listener weights have unit squared norm and every pan pair has unit
power, the decoder's Frobenius norm is one and its spectral norm cannot exceed
one. No fitted gain correction is required for passivity.

A listener near a side wall may receive a real lateral energy bias because that
wall has greater coupling. A centred listener in a left/right symmetric room
must have exactly mirrored decode coefficients.

### What the late tail should preserve

The late field is not expected to keep a source hard-panned indefinitely.
Instead:

- the first late energy inherits bias from the source's six wall sends;
- repeated FDN mixing reduces that source-specific bias;
- mature tail energy is balanced for a symmetric room and centred listener;
- the channels remain decorrelated enough to create width and envelopment;
- listener position continues to affect wall delays and wall energy;
- the response remains deterministic for a static scene.

Mirror symmetry for ER is a sample-for-sample requirement. For the mixed late
tail, mirrored scenes are required to match energy, decay, correlation, and
lateral-bias statistics; their individual samples need not be identical.

### Contingency for insufficient late decorrelation

The six-wall decoder is the first production implementation. If its measured
late interchannel correlation fails the limits below, add one explicitly
diffuse residual stereo pair from the ten-dimensional complement of
`WallProjection`.

Such residual rows must:

- be orthogonal to the six wall columns and to each other;
- have equal energy and deterministic coefficients;
- be described as geometry-independent diffuse energy, not virtual directions;
- be mixed at a fixed, documented gain;
- leave the combined output operator passive after normalization.

Do not add this residual unless the six-wall implementation fails an objective
correlation or listening requirement.

## Width control

The physical renderer is defined before Width. Width remains an intentional
speaker-output transformation on the complete wet signal, preserving the
current product behaviour:

```text
mid  = (L + R) / sqrt(2)
side = width * (L - R) / sqrt(2)

Lout = (mid + side) / sqrt(2)
Rout = (mid - side) / sqrt(2)
```

The existing mapping from the normalized control to 0-150% width remains:

```text
width = 1.5 * smoothstep(control)
```

The factory control value maps to unity width. Zero Width makes the complete
wet signal mono; maximum Width exaggerates the base room image. Spatial
correctness tests must measure the base renderer at unity width, before this
creative transformation.

## Diffusion control

### Required semantic change

Diffusion must control both:

1. how strongly each Hadamard butterfly mixes its pair; and
2. the temporal span of both velvet-delay stages.

Keep the existing delay-span law:

```text
u         = smoothstep(diffusion)
spanScale = 0.20 + 1.30 * u
```

multiplied by the current room scale.

Replace the fixed 45-degree butterfly with

```text
[ y0 ]   [ cos(theta)   sin(theta) ] [ x0 ]
[ y1 ] = [ sin(theta)  -cos(theta) ] [ x1 ]
```

and initially use

```text
thetaMin = pi / 16
theta    = thetaMin + u * (pi / 4 - thetaMin)
```

At minimum the network is mildly mixed rather than reduced to independent comb
loops. At maximum it is the present normalized Hadamard butterfly.

Use the same `theta` law for every butterfly layer and both coefficient
flavours for the first implementation. If echo-density tests reveal a
non-monotonic region, change the per-layer angle schedule explicitly and record
it as part of the topology; do not repair it with fitted output EQ or gain.

### Energy and automation rules

The parameterized pair matrix is orthogonal for every `theta`. Interpolate
`theta` itself during automation, so every intermediate matrix remains
orthogonal. Do not linearly crossfade the outputs of two diffusion matrices.

Signed permutations remain passive routing. They do not count as mixing when
describing the Diffusion control.

Orthogonality guarantees that the instantaneous matrix does not change vector
energy. It does not by itself guarantee unchanged measured RT60 once unequal
delays and line-dependent loss filters are included. RT60 invariance must
therefore be measured, not inferred.

Velvet tap changes continue to crossfade fixed integer taps without clearing
buffers. No Diffusion automation path may call `Reset()`.

## Distance and level policy

After this correction, automatic distance behaviour consists only of:

- bounded, room-scale-normalized direct gain;
- physical inverse-distance spreading on every ER image path;
- existing wet Size calibration;
- source-specific tail send derived from ER handoff power;
- normalized source and listener wall direction vectors.

Remove `PositionBalanceGains`, `currentEarlyPositionGain_`,
`currentTailPositionGain_`, and their targets and smoothing state.

`EARLY` and `TAIL` remain exact user decibel trims around this baseline. They
must not alter direct sound, network feedback poles, or decay time.

If listening tests later demonstrate that the direct-to-reverberant distance
range remains inadequate, add only a bounded documented residual after
measuring the terms above. Do not restore an equal-and-opposite ER/tail proxy by
default.

## Required verification

All rendered tests use the production 48 kHz sample rate and the complete
16-line topology.

Use the following definitions consistently:

```text
EL = sum(L[n]^2) over the measurement window
ER = sum(R[n]^2) over the measurement window

lateralBias = (EL - ER) / max(EL + ER, epsilon)

correlation = sum(L[n] * R[n]) /
              sqrt(max(sum(L[n]^2) * sum(R[n]^2), epsilon))

effectiveLineCount = (sum(lineEnergy[i]))^2 /
                     sum(lineEnergy[i]^2)
```

Apply the specified band-pass filters before calculating band correlation or
band energy. Window positions are measured from the logical direct reference,
not from the host input sample.

### Direct tests

- Hard-left, centred, and hard-right test vectors match the analytical pan law.
- A source left of the listener produces greater direct energy on the left.
- Mirroring source and listener X positions swaps direct L/R samples exactly.
- Moving source and listener together while preserving their relative vector
  preserves direct pan and normalized distance gain.
- Uniform Size scaling at fixed normalized positions changes direct gain by
  less than 0.25 dB.
- The bounded distance law never exceeds +6 dB or falls below -12 dB.
- Multiple sources are independently panned before summation.
- `CentredLegacy` produces identical L/R dry samples.
- Direct output begins on the input sample with no implementation delay.

### ER tests

- Analytical first-order paths match measured excess-delay tap times within one
  fractional-delay interpolation tolerance.
- ER output begins at `round(tExcess * sampleRate)` samples before user Pre-delay
  is added.
- Every isolated path preserves power through its pan.
- A centred path is exactly channel-balanced.
- A left path has positive left-energy bias; a right path has the opposite sign.
- Mirroring the complete room in X swaps ER kernels within numerical tolerance.
- Front/rear moves alter the FIR only through calculated path geometry and
  material encounters.
- The ER envelope is monotonic over the handoff interval and identical in both
  channels.
- Geometry updates crossfade prepared FIR banks without clearing convolution or
  input history.

### ER/late transition tests

- Late onset occurs no earlier than 2 ms before `handoff.finalStartSeconds` and
  no later than 2 ms after it.
- Windowed total wet energy across the handoff contains no dip or step greater
  than 1.5 dB relative to adjacent equal-length windows.
- Each four-band handoff-energy step is below 2 dB.
- EARLY and TAIL controls remain exact decibel offsets around the automatic
  baseline.
- Source alignment-delay automation preserves already stored late energy.

### Late stereo tests

Measure windows relative to the direct reference and exclude Width processing.

- A symmetric room and centred listener have L/R energy imbalance below 0.5 dB
  after 250 ms.
- In 125-500 Hz, absolute normalized interchannel correlation is below 0.80
  after 250 ms.
- In 500 Hz-4 kHz, absolute normalized interchannel correlation is below 0.50
  after 250 ms.
- Mono-summed late energy must not fall more than 6 dB below the average stereo
  channel energy in any test band.
- For an off-centre source, initial late lateral bias has the same sign as the
  source. After `max(250 ms, 2 * handoff.finalEndSeconds)`, its magnitude is at
  most 25% of the initial bias in a symmetric room.
- Mirrored scenes match late L/R energy, decay time, and band correlation after
  channel swapping within explicit numerical tolerances.
- Moving the listener changes wall timing and energy without changing feedback
  poles or RT60.
- Zero Width is exactly mono, unity Width equals the base decode, and Width
  changes side energy monotonically without changing mid energy.

Correlation limits are initial production gates. If listening validation shows
that one is unnecessarily strict or loose, change it once with recorded impulse
responses and rationale; do not omit the assertion.

### Diffusion tests

- The lossless VFM conserves energy at diffusion 0, 0.25, 0.5, 0.75, and 1.
- Effective line count over the first two network recurrences rises
  monotonically as Diffusion increases.
- Normalized echo density and temporal support rise monotonically at the same
  control points.
- Mid-band RT60 changes by less than 5% across the Diffusion sweep.
- Low- and high-band RT60 change by less than 7% across the sweep.
- Modal prominence and low-frequency periodicity remain within the existing
  limits over the complete Size, Decay, Damping, and Diffusion grid.
- Maximum Size and Decay remain finite and decaying at every Diffusion value.
- Diffusion automation retains all delay buffers and creates no click, large
  energy hole, or Doppler sweep.

## Code change map

| File | Required change |
| --- | --- |
| `src/tfdsp/early_reflections.hpp` | Use horizontal lateral direction; place FIR taps in excess time; remove latency-compensation response fields; split prepared FIRs into zero-latency heads and partitioned tails. |
| `src/tfdsp/early_reflections_worker.*` | Publish the revised FIR and handoff representation without doing new work on the audio thread. |
| `src/tfdsp/room_reverb.hpp` | Produce stereo direct and wet frames; add direct mode and direct-distance law; remove `PositionBalanceGains`; add per-source handoff-alignment delays. |
| `src/tfdsp/late_reverb.hpp` | Retain six wall ports; replace `SidePattern` with the physical wall pan decoder; expose or calculate intrinsic late onset needed for alignment. |
| `src/tfdsp/velvet_feedback_matrix.hpp` | Pass smoothed `theta` through every butterfly layer while retaining the existing velvet-span transition. |
| `src/tfdsp/reverb_output.hpp` | Accept a stereo direct frame in both output-mix overloads. |
| `src/TfReverb.cpp` | Remove scalar dry summation; consume `RoomReverbFrame`; serialize and display `DirectMode`; update parameter descriptions. |
| `tests/early_reflections_tests.cpp` | Add analytical excess-time, horizontal-pan, mirror, and zero-latency head/tail tests. |
| `tests/late_reverb_tests.cpp` | Add physical decoder, onset, band-correlation, lateral-bias decay, and variable-mixing tests. |
| `tests/room_reverb_tests.cpp` | Add complete direct/ER/tail timing, positioned polyphony, distance, transition, Width, and migration tests. |

## Implementation sequence

Implement and review the correction in the following order:

1. Add failing analytical pan, mirror, excess-time, correlation, lateral-bias,
   and diffusion-strength tests.
2. Introduce `RoomReverbFrame`, stereo direct mixing, and `DirectMode` while
   retaining a zero-latency direct path.
3. Convert ER FIR placement and handoff analysis completely to excess time;
   implement the zero-latency FIR head/tail split and remove convolution-latency
   compensation state.
4. Remove `PositionBalanceGains` and position-gain smoothing. Retain the
   handoff-derived scalar late send.
5. Add source-specific late alignment delays and verify ER/late overlap across
   the Size and Aspect ranges.
6. Replace `SidePattern` with the normalized six-wall constant-power decoder.
7. Measure late correlation before deciding whether the optional diffuse
   complementary output pair is needed.
8. Parameterize the VFM butterfly angle and add state-preserving smoothing.
9. Update module serialization, menu text, parameter descriptions, current
   design documentation, and patch migration.
10. Run native tests, automation tests, benchmarks, and listening comparisons at
    48 kHz before changing any FDN coefficient set.

Once these stages pass, the corrected renderer becomes the sole production
spatial path. `CentredLegacy` affects only direct-sound presentation; it must not
retain the old ER timing, side pattern, distance proxy, or Diffusion behaviour.

## Comparison with representative open-source reverbs

The design deliberately differs from common musical stereo reverbs:

- Freeverb creates width by feeding a mono sum into two networks with different
  delay lengths, then crossmixing their outputs. It does not preserve room or
  source direction.
- Dragonfly/Freeverb3 keeps stereo tank halves, cross-feedback, distinct output
  taps, and variable allpass diffusion. Its spatial image is designed rather
  than geometric.
- CloudSeed uses separate randomized left/right engines and controls how similar
  their seeds and inputs are. Its late width is decorrelation, not a room decode.
- zita-rev1 uses a fixed stereo input/output matrix around an eight-line FDN and
  supplies a separately designed four-channel Ambisonic output mode.
- Steam Audio preserves a directional Ambisonic reflection field until a final
  speaker or binaural decoder.

TfReverb sits between the musical and physical approaches: it has real
image-source ER and meaningful six-wall late connections, but only two-speaker
output. Keeping the six wall ports and decoding their true lateral directions
preserves the geometry the topology actually owns. Treating the mixed FDN lines
as 16 virtual directions would not.

## References

- GRAME, `faustlibraries/reverbs.lib`: Freeverb, Dattorro, and zita-rev1 stereo
  and Ambisonic implementations.
  <https://github.com/grame-cncm/faustlibraries/blob/master/reverbs.lib>
- Michael Willis and Freeverb3, Dragonfly Reverb source.
  <https://github.com/michaelwillis/dragonfly-reverb>
- Ghost Note Audio, CloudSeed Core source.
  <https://github.com/GhostNoteAudio/CloudSeedCore>
- Valve, Steam Audio Reflection Effect. Reflections remain Ambisonic before
  speaker or binaural decoding.
  <https://valvesoftware.github.io/steam-audio/doc/capi/reflections-effect.html>
- C. Kirsch, J. Poppitz, T. Wendt, S. van de Par, and S. D. Ewert, "Spatial
  Resolution of Late Reverberation in Virtual Acoustic Environments," *Trends
  in Hearing*, 25, 2021. Its 12-24 virtual-source result concerns spatial
  reproduction over a spherical loudspeaker array, not ordinary stereo output.
  <https://doi.org/10.1177/23312165211054924>
- E. De Sena, H. Hacihabiboglu, Z. Cvetkovic, and J. O. Smith, "Efficient
  Synthesis of Room Acoustics via Scattering Delay Networks," *IEEE/ACM TASLP*,
  23(9), 2015. This is the reference model for geometry-connected passive wall,
  source, and receiver nodes.
  <https://doi.org/10.1109/TASLP.2015.2438547>
