# TfReverb current design

TfReverb is a mono/polyphonic-input, stereo-output rectangular-room reverb. It
combines deterministic geometric early reflections with one passive late-field
architecture. There are no alternate, decay-dependent, or legacy late
topologies in the production path.

## Signal path

```text
mono sources with XYZ positions
        |
shared wet pre-delay
        |
        +-- rectangular image-source early reflections --+
        |                                                |
        +-- six position-dependent wall sends            |
             -> 16 geometry-scaled main delays           |
             -> three-band loss per delay segment        |
             -> four-coordinate filtered octave blend    |
             -> geometry-scaled two-stage velvet matrix -+
                       ^                                  |
                       +-- optional delay modulation -----+
                                                          |
geometry- and distance-matched ER/tail balance -> stereo decode -> wet filters -> width
```

The early-reflection worker builds four-band FIRs away from the audio thread.
Prepared convolution banks crossfade on the audio thread. Source position,
listener position, room dimensions, and material damping affect the geometric
response. Diffusion does not change the number or truncation horizon of image
sources. Size derives length, width, and height; Aspect reshapes the floor while
preserving its area. A characteristic-length gain compensates inverse-distance
spreading so Size does not also become an unintended wet/dry control.

## Late field

The late field has 16 main delays. Their selected relative ratios are
normalized to unit mean, so changing coefficient flavour cannot redefine the
Size control. Their absolute mean delay is derived from the room mean free time

```text
t_mean = 4 V / (c S)
```

where `V` is room volume, `S` is surface area, and `c` is the speed of sound.
This keeps the early-reflection geometry and recursive late-field time scale
coherent. The room mean-free-time ratio also scales both velvet delay banks, so
no fixed recursive delay span can dominate the Size control. Main-delay changes
use dual read heads and a convex crossfade; velvet stages crossfade fixed integer
taps. Control movement does not sweep a read head through stored audio.

The feedback operator is a two-stage velvet feedback matrix

```text
U2 D2(z) U1 D1(z) U0
```

Each `U` is a fixed normalized Hadamard butterfly followed by a signed
permutation. Each `D` is a diagonal bank of pure delays. At a static control
setting the operator is paraunitary. Diffusion scales the temporal span of both
`D` banks from a nonzero safe minimum to 1.5 times their room-scaled nominal
span; it does not weaken or rotate the feedback transforms. Runtime changes
crossfade between integer delay taps and remain passive.

The deterministic heuristics in `src/tfdsp/late_reverb_coefficients.hpp`
define the immutable Base flavour. A separately exported fit defines the
Optimized flavour without changing the topology or control laws. The module's
**Late-tail FDN** context menu selects **Base FDN** or **Optimized FDN**;
Optimized FDN is the default and Base FDN remains available as the fixed
comparison reference.
Switching crossfades the main-delay reads, velvet taps, and passive transform
blend while retaining every delay buffer, so the same live tail can be used
for a direct comparison.

The accepted Optimized flavour deliberately retains Base's velvet delays,
signed permutations, and signs. An unconstrained discrete candidate improved
the aggregate frequency-domain objective but failed the native maximum-room,
minimum-Diffusion 20--80 Hz periodicity limit. The production fit therefore
changes only the normalized main-delay ratios, each within 0.005 of its Base
heuristic, and uses the midpoint between Base and the converged fit. This makes
the context-menu comparison specifically about late-mode spacing rather than
two different scattering networks.

## Decay and damping

Decay is the mid-band RT60 and ranges from 0.25 to 8 seconds. Every physical
delay segment, including both velvet stages, has a complementary three-band
attenuation filter:

- low: below approximately 220 Hz;
- mid: approximately 220 Hz to 3.5 kHz;
- high: above approximately 3.5 kHz.

At zero Damping all bands use the selected Decay. At maximum Damping, low-band
RT60 is 72% and high-band RT60 is 20% of the mid-band value. Equal band gains
reconstruct the input exactly, so this block controls decay rather than adding
an unrelated static equalizer.

## Diffusion, modulation, and shimmer

Diffusion controls short-term temporal scattering. The transforms remain fully
mixed at every setting, while the two velvet delay spans change substantially.
Minimum Diffusion therefore remains stable and avoids exposed comb filters;
maximum Diffusion produces a much wider cloud around every recurrence.

Modulation applies independent, slow random motion of the 16 main-delay read
positions. The lower 35% of the knob is a deliberately static region. Above it,
depth rises cubically to 0.25 ms peak excursion.

Shimmer is part of the recursive late-room loop, as it is in the established
feedback-shimmer topology. It projects the attenuated 16-line late field onto
four orthonormal Walsh coordinates. Each is high-pass filtered, shifted up one
octave by an independently seeded eight-grain shifter, and darkened by two
low-pass stages at `min(6.5 kHz, 0.2 * sample rate)`. The knob adds a bounded
shifted return of up to 85% within those four coordinates without replacing
the original late field. The other 12 coordinates are untouched.

The reconstructed 16-line field then passes through the complete two-stage
velvet feedback matrix and the main room delays before it can reach the stereo
decoder. There is no direct pitch-shifter output, separate shimmer tank, or
decoder-side shifted layer. Consequently every octave return is scattered by
the production diffuser, while orthonormal projection confines the bounded
return to a known subspace rather than creating an independent tank. Grain
durations stay equal so their uniformly staggered Hann windows maintain
constant overlap; only look-back and pitch receive small deterministic random
offsets.

## Position and stereo

Every source is connected to six wall ports using distance delays and
solid-angle-derived normalized weights. The listener uses an equivalent set of
six output connections. The stereo decoder is factorized into orthogonal mid
and side wall modes, guaranteeing equal left/right decoder energy without a
fitted balance correction.

The early-reflection handoff analysis selects a source-dependent baseline late
send for a continuous spectral transition. A separate source-listener distance
law preserves the perceptual distinction between the discrete early field and
the approximately diffuse late field: nearby pairs favor ER definition, while
distant pairs favor the tail. The law uses physical metres, follows a bounded
20-log distance ratio around the factory placement, and applies equal-and-
opposite gains so it primarily changes balance rather than raw wet level.
EARLY and TAIL remain user trims around that automatic baseline. STEREO WIDTH
is an output mid/side operation and does not change room geometry.

Wet pre-delay, room-size delays, velvet-stage delays, and source/listener wall
delays all retain their buffers during automation. Time changes use old/new
read heads and smooth crossfades; no parameter path calls `Reset()`.

## Required invariants

The native test suite verifies that:

- all rendered-audio regressions run at the 48 kHz production sample rate;
- each static VFM setting conserves energy in both coefficient flavours;
- Diffusion decisively expands temporal scattering support;
- Size monotonically increases every room dimension and the complete recursive
  delay span;
- Size compensation prevents a large wet-level inversion;
- Damping shortens high-frequency decay much more than low-mid decay;
- maximum Size and Decay remain finite, decaying, dense, and non-periodic in
  both coefficient flavours;
- automated size changes crossfade without Doppler sidebands;
- Shimmer enters only through the recursive late-field path, increases octave
  energy monotonically, remains spectrally diffuse and free of grain-boundary
  spikes, reaches a useful octave/fundamental ratio, and cannot destabilize the
  room tail;
- source-listener distance monotonically moves the automatic ER/tail balance;
- control automation and coefficient-flavour changes preserve stored tail and
  pre-delay history;
- the default stereo response has no left/right energy bias.

The coefficient search reproduces this exact architecture and must retain these
invariants. It screens delay ratios, delay stratification, and signed
permutations, but must not introduce a second runtime topology or replace the
deterministic control semantics.
