# TfReverb spatial and diffusion correction

## Status and scope

This document specifies the correction target for TfReverb's source/listener
positioning, stereo early reflections, late spatial field, ER/late transition,
and Diffusion control. It records intended production behaviour; the changes
described here are not yet implemented.

The target remains a mono or polyphonic-source, stereo-output room reverb for
ordinary speaker and production use. Binaural rendering, HRTFs, head tracking,
and literal behind-the-listener localization are out of scope.

Front/rear differences must come from the rectangular room geometry: different
path lengths, wall encounters, attenuation, and reflection patterns. The engine
must not add a synthetic rear EQ, pinna model, or other front/back heuristic.
The room-plan top edge defines listener-forward only so left/right directions
and geometry have a stable coordinate frame.

## Why correction is required

The geometry engine currently calculates considerably more spatial information
than the stereo signal path uses:

- early-reflection paths retain a three-dimensional arrival direction, but the
  renderer reduces it to one frequency-independent equal-power pan;
- the direct signals of every source are summed to centred mono before the
  dry/wet mix, so source position cannot localize the dominant direct sound;
- ER taps are placed at absolute propagation time while the direct signal is
  emitted without its corresponding propagation delay;
- late source/listener connections remove common wall delay and normalize wall
  energy, then a fixed equal-energy side pattern replaces a directional stereo
  decode;
- the distance law changes ER/tail balance without acting on the direct sound;
- Diffusion leaves every feedback transform fully mixed and changes only the
  internal velvet-delay span.

The existing tests mostly prove that position changes some samples and that
the output is not exactly mono. They do not establish correct lateral sign,
mirror symmetry, direct-relative timing, directional persistence, stereo
coherence, or perceptually meaningful diffusion behaviour.

These are architectural and verification problems, not coefficient-tuning
problems. No optimization should be used to compensate for them.

## Target signal path

```text
independent mono sources with XYZ positions
        |
        +-- positioned direct renderer -----------------> stereo direct
        |
        +-- shared user pre-delay
             |
             +-- direct-relative image-source ER
             |       -> directional stereo decode ------+
             |                                           |
             +-- six geometry-dependent wall sends      |
                    -> 16 room-scaled main delays        |
                    -> passive multiband loss            |
                    -> variable-strength two-stage VFM   |
                    -> 12--16 directional late outputs   |
                    -> stereo speaker decode ------------+
                                                         |
                         automatic ER/tail balance -> wet filters
                                                         |
stereo direct ----------------> direct/wet mix <----------+
                                      |
                                 output level
```

Direct sound, discrete reflections, and diffuse tail have different spatial
roles:

- direct sound is strongly positioned;
- the first and strongest ERs remain strongly directional;
- higher-order ERs lose directional specificity while approaching the mixing
  interval;
- the beginning of the tail may retain a weak source-dependent bias;
- the established tail is spacious and mostly source-independent.

It is correct for a diffuse late field to lose the direction of individual
reflections. It is not correct to lose its lateral energy distribution and
interchannel coherence. The late renderer therefore preserves spatial
statistics rather than a permanently panned source direction.

## Coordinate and orientation convention

Source and listener positions remain normalized XYZ coordinates inside the
room. In the two-dimensional room plan:

- increasing X moves right;
- the top edge is listener-forward;
- increasing/decreasing Y changes front/rear geometry according to the panel's
  established orientation;
- Z continues to affect physical path lengths and wall interactions even though
  the main panel exposes only the plan view.

All lateral decisions must use source or path direction relative to the
listener, not absolute source X. Moving the source and listener together while
preserving their relative geometry must preserve the direct stereo placement.

No additional front/back coloration is applied. A source moved from front to
rear may sound different because the set of image paths and wall encounters is
different, but the stereo output does not claim to place sound behind the
listener.

## Direct path

The direct path must remain separate per source until stereo decoding. The
current polyphonic-to-mono sum before the output mix must be removed from the
positioned mode.

For each active source:

1. calculate the source-to-listener vector in metres;
2. derive its normalized lateral component or listener-relative azimuth;
3. apply a constant-power stereo pan;
4. sum the resulting stereo source pairs.

Artificial interchannel delay is not added in speaker-stereo mode. Constant-
power amplitude panning is predictable, mono-compatible, and does not create
comb filtering when the two channels are summed.

The recommended product behaviour is **Position direct sound**, in which the
dry side of the Mix control is the positioned stereo direct signal. A context
option may retain **Centred direct sound** for conventional insert use and
backward compatibility. This choice must be explicit; the engine must not call
a centred mono sum physically positioned sound.

Distance gain on the direct path requires a bounded reference-distance law.
It must avoid singular gain when source and listener coincide and must not make
moving the room pad an uncontrolled output-level gesture. The final law should
be calibrated together with ER spreading and the existing overall output
level, rather than layered on top of the current ER/tail proxy unchanged.

## Early-reflection correction

The rectangular image-source enumeration, material encounter counts,
frequency-band absorption, air loss, and worker-thread FIR construction remain
the basis of the ER engine.

### Timing reference

The positioned direct sound is the zero-time reference. Each reflection tap is
placed from its excess path length:

```text
t_ER = (d_reflection - d_direct) / c
```

The engine already calculates this value as `excessDelaySeconds`. Absolute
propagation time must not be used for a reflection while its corresponding
direct sound remains at zero latency. The user Pre-delay is applied afterward
as an explicit creative offset to the wet paths.

### Stereo direction

Every retained path is decoded from its physical listener-relative arrival
direction. For stereo speakers, only the lateral projection controls channel
balance. Front/rear and elevation continue to influence geometry but do not
receive invented localization processing.

Constant-power path gains must satisfy equal power before distance and material
losses. Left/right room mirror cases must produce mirrored kernels. A centred,
geometrically symmetric path must be exactly balanced.

### Directional transition with reflection order

The first and strongest reflections carry the clearest directional
information. Directional specificity should decrease smoothly as the response
approaches the detected mixing interval. This can be implemented by blending
each path's physical stereo pan toward an energy-normalized diffuse stereo
distribution as a function of reflection order and/or normalized time to the
mixing interval.

The blend must not alter total path power. It must also be deterministic, so a
static scene produces a static FIR and geometry updates can crossfade without
image wander.

## Distance and ER/tail balance

Once direct and ER timing and gains share one reference, the existing
`PositionBalanceGains` law must be reassessed rather than preserved by default.
It currently applies equal-and-opposite gains only to ER and tail and therefore
does not implement a genuine direct-to-reverberant distance cue.

The corrected balance should arise primarily from:

- bounded direct-distance gain;
- image-path distance spreading;
- source-dependent ER energy at the mixing interval;
- geometry-dependent late injection.

EARLY and TAIL remain user trims around the automatic baseline. Any remaining
perceptual distance correction must be a documented bounded residual, measured
after the physically connected terms above, rather than a second full-strength
distance law.

## ER/late transition

The detected ER mixing interval must control more than ER truncation and a
scalar late-send gain.

The production transition must:

- taper directional ER energy over the detected start/end interval;
- arrange the FDN buildup so dense late energy arrives through the same
  interval;
- retain overlap so there is neither an energy hole nor a hard boundary;
- preserve the spectrum measured at the handoff;
- avoid an early diffuse wash that masks the first directional reflections.

Possible implementations include delaying or pre-seeding each source's late
injection so the first useful FDN output meets the mixing interval, or applying
an energy-preserving buildup envelope to the late output. The choice must be
made from measured impulse responses, not solely from implementation
convenience.

## Late spatial field

### What should and should not survive diffusion

The mature tail should not remain hard-panned to the source. Diffusion is
expected to destroy individual path direction and make the tail progressively
source-independent.

The following must survive as controlled stereo properties:

- balanced total left/right energy for a symmetric room;
- useful interchannel decorrelation;
- lateral energy and envelopment;
- deterministic response for a static scene;
- room/listener-dependent anisotropy where geometry or wall energy supports it;
- a gradual loss of initial source-direction bias rather than an instantaneous
  arbitrary remap.

### Directional virtual outputs

The 16 FDN lines should feed 12--16 decorrelated virtual reverb outputs with
fixed room-relative directions. Using all 16 is the preferred starting point
because the signals already exist and avoids reducing the internal spatial
resolution before stereo decoding.

The virtual-output projection must be passive or explicitly normalized. Its
direction set should cover the room symmetrically. Horizontal directions are
constant-power panned into stereo; vertical directions contribute centrally or
with symmetric lateral weighting. Front and rear directions use the same
speaker-stereo lateral law and differ only through the room field that reaches
them.

Listener position and directional wall energy may weight the virtual outputs,
but must not introduce an unexplained fixed `SidePattern`. In a uniform room,
the mature late field should approach a balanced diffuse stereo field even when
the source begins off-centre.

The STEREO WIDTH control remains a speaker-output control. It may adjust the
late field's side energy after a physically meaningful base decode, but it must
not substitute for that decode or change the underlying room geometry.

## Diffusion correction

### Current problem

The current VFM applies fully mixing Hadamard transforms at every Diffusion
setting. The control changes only the two internal velvet-delay spans. This is
stable, but it does not control scattering strength and explains why low and
high settings can sound insufficiently distinct.

### Target operator

Diffusion controls both orthogonal mixing strength and temporal scattering
span. Each butterfly becomes a parameterized energy-preserving rotation. A
representative pair operation is

```text
[ y0 ]   [ cos(theta)   sin(theta) ] [ x0 ]
[ y1 ] = [ sin(theta)  -cos(theta) ] [ x1 ]
```

with `theta` increasing smoothly from a safe mild-scattering value to the full
design angle. Signed permutations remain passive routing operations; they must
not be counted as signal mixing merely because they move a line index.

The control law should have these endpoints:

- minimum: mild but nonzero scattering and short velvet span, retaining
  separated texture without exposing unstable or strongly periodic comb loops;
- default: a dense, natural room buildup;
- maximum: full scattering and the maximum validated velvet span.

Every static operator remains orthogonal/paraunitary before its explicit loss
filters. Diffusion therefore cannot change loop energy or RT60 by accident.
Runtime changes use state-preserving coefficient/tap transitions and never
clear delay buffers.

The Base and Optimized FDN flavours retain the same Diffusion law. Their menu
comparison remains a coefficient-flavour comparison, not a topology change.

## Objective verification

All rendered-audio tests run at the 48 kHz production sample rate. No
downsampled network or shortened physical control law may stand in for the
production engine.

### Direct and ER invariants

- a left source produces greater early/direct energy on the left;
- a right source is the numerical mirror of the corresponding left source;
- moving the listener across a fixed source reverses lateral placement;
- a centred source in a symmetric scene remains channel-balanced;
- moving source and listener together while preserving relative geometry
  preserves direct lateral placement;
- analytical first-order image paths match measured excess-delay tap times;
- reflection panning preserves path power;
- front/rear moves change the IR only through their calculated geometry;
- multiple sources remain independently positioned through the direct and ER
  paths;
- geometry updates crossfade without clearing input or convolution history.

### ER/late invariants

- ER directionality decreases toward the mixing interval;
- tail buildup overlaps the ER fade without a statistically significant energy
  hole or step;
- the handoff spectrum is continuous within explicit band tolerances;
- source distance produces a monotonic, bounded direct/ER/tail relationship;
- EARLY and TAIL trims remain exact decibel offsets around the automatic
  baseline.

### Late stereo invariants

- virtual-output energy is bounded and normalized;
- left/right energy is equal for symmetric rooms and mirrored for mirrored
  scenes;
- pairwise virtual-output correlation and final stereo interchannel
  correlation remain within explicit limits over several time windows and
  frequency bands;
- source-direction bias decreases with time while late spatial extent remains;
- Width changes side energy monotonically without changing RT60 or room
  geometry.

### Diffusion invariants

- the lossless VFM conserves energy at minimum, default, maximum, and
  intermediate settings;
- echo density and temporal support rise monotonically with Diffusion;
- Diffusion produces a decisive response difference at several control points;
- RT60 and multiband damping remain invariant within tolerance;
- modal prominence and low-frequency periodicity remain bounded across the
  complete Size, Decay, Damping, and Diffusion grid;
- static controls remain finite and decaying at maximum Size and Decay;
- automation retains every delay buffer and produces no large amplitude hole,
  click, or Doppler sweep.

Tests that merely assert that two responses differ are insufficient for spatial
or diffusion controls. Each test must assert the sign, monotonic direction,
symmetry, physical timing, correlation, density, or stability property that the
control is intended to provide.

## Implementation sequence

1. Add directional metrics and failing regression tests for the current direct,
   ER, late-stereo, timing, and Diffusion behaviour.
2. Separate and stereo-position every direct source; define the explicit
   positioned/centred direct policy.
3. Change ER tap placement to direct-relative excess delay and implement
   physically directed stereo kernels.
4. Replace the distance proxy with a calibrated direct/ER/tail relationship and
   align FDN buildup with the detected mixing interval.
5. Replace the fixed two-row late `SidePattern` decoder with 12--16 normalized
   directional virtual outputs and a speaker-stereo decode.
6. Parameterize the VFM butterflies so Diffusion controls real scattering
   strength as well as velvet span.
7. Run the complete native 48 kHz response, stability, automation, and
   performance suites before listening comparison.
8. Fine-tune only after the architecture and invariants pass. Optimization may
   improve modal flatness or correlation, but must not define basic spatial or
   control semantics.

Once the corrected design passes objective and listening validation, it becomes
the single production path. The present spatial decoder and delay-reference
shortcuts should be removed rather than retained as an undocumented alternate
architecture.

## References

- T. Wendt, S. van de Par, and S. D. Ewert, “A Computationally-Efficient and
  Perceptually-Plausible Algorithm for Binaural Room Impulse Response
  Simulation,” *Journal of the Audio Engineering Society*, 62(11), 2014.
  The relevant architectural point is the hybrid image-source/FDN split and the
  use of spatially distributed late outputs; this specification does not adopt
  its binaural renderer.
  <https://doi.org/10.17743/jaes.2014.0042>
- C. Kirsch, J. Poppitz, T. Wendt, S. van de Par, and S. D. Ewert, “Spatial
  Resolution of Late Reverberation in Virtual Acoustic Environments,” *Trends
  in Hearing*, 25, 2021. The study found that 12 spatially distributed virtual
  late sources could approximate higher resolutions in its tested isotropic
  conditions, motivating retention of multiple directional late outputs before
  final stereo decode.
  <https://doi.org/10.1177/23312165211054924>
- E. De Sena, H. Hacıhabiboğlu, Z. Cvetković, and J. O. Smith, “Efficient
  Synthesis of Room Acoustics via Scattering Delay Networks,” *IEEE/ACM
  Transactions on Audio, Speech, and Language Processing*, 23(9), 2015. This
  provides a reference for geometry-connected passive scattering and
  directional source/receiver treatment.
  <https://doi.org/10.1109/TASLP.2015.2438547>
- Valve, “Steam Audio Reflection Effect.” Its reflected field is retained in an
  intermediate directional representation before output decoding, illustrating
  the distinction between losing individual late rays and discarding the late
  field's spatial distribution.
  <https://valvesoftware.github.io/steam-audio/doc/capi/reflections-effect.html>
