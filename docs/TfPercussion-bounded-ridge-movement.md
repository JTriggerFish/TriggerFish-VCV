# Bounded ridge movement — prototype

The purpose is moving, irregular metallic detail without necessarily replacing
each coherent oscillator with a broad noise band. This is a constructive
synthesis experiment, not an identified physical model of this reference gong.
Existing phase blur and slow pitch wander remain independently available for
comparison. At the user's request, the gong preset now loads the recommended
blend: movement 1 rad, speed 100 changes/s, packet sharing 0.25, together with
the existing phase blur 0.035 ERB and blur tilt -0.5. Only those three movement
values changed; depth zero restores the previous sound. Other presets stay off.

## Controls and integration

The workbench exposes **Modal body → Ridge movement · experiment**:

| UI / JSON | Meaning | Default |
|---|---|---|
| Ridge movement / `field_motion_depth` | Maximum bounded phase displacement, 0–3 radians | 0: bypass |
| Movement speed / `field_motion_rate` | Mean random trajectory knot rate, 0.1–200 changes/s | 40 |
| Packet sharing / `field_motion_sharing` | Mix of independent mode movement and shared packet movement | 0.5 |

All three values pass from JSON through the existing WASM descriptors to
`CrashCymbalFitParameters::fieldMotion`, then to the prepared modal field.
They are visible and saved. Old complete snapshots gain the off-by-default
triple explicitly on import; partial triples are rejected. Crash, ride, gong
and hi-hat factory snapshots store the triple explicitly. No new build
dependency is introduced for Rack; the web workbench remains optional.

The gong preset loads depth **1**, speed **100**, sharing **0.25**. Returning
depth to zero restores the previous renderer path. This is an audition starting
point, not a claimed reference calibration. Reducing Phase blur helps isolate
the distinction, but removing all blur also removes useful diffuse energy.

## Algorithm

Each mode and each packet has a bounded random trajectory between independent
uniform targets in $[-1,1]$. Quintic smoothstep interpolates each segment:

$$s(t)=6t^5-15t^4+10t^3,\quad 0\leq t\leq1.$$

Initial segment positions differ; rates vary uniformly between 0.75 and 1.25
times the speed control at each new target. There is no common periodic LFO.
For mode $i$ in packet $p$, with sharing $c$ and depth $d_i$:

$$q_i[n]=d_i\big((1-c)u_i[n]+c\,v_p[n]\big).$$

This convex mix stays within the stated maximum depth. Its variance is smaller
in the middle of the sharing range than at either end; the tooltip states that
the middle is gentler. It is not an RMS-normalized crossfade or a correlation
coefficient. Different packets have independent shared trajectories.

The extra oscillator rotation is the **difference** of this displacement:

$$z_i[n+1]=r_i\exp\left(j\left[\omega_i+
q_i[n+1]-q_i[n]\right]\right)z_i[n]+b_i x[n].$$

Thus the added phase telescopes to a bounded difference, rather than walking
away indefinitely. In an isolated unforced oscillator, modest depth retains a
coherent carrier plus moving sidebands. Rotation preserves the complex-state
norm algebraically; the existing damping radius is unchanged. Input can still
interfere constructively or destructively with stored vibration, and the
existing energy cascade still operates. This does not guarantee unchanged
observed loudness, an identical fitted envelope, or no spectral sidebands.

The implementation uses Cayley rotations and a small-angle tangent expansion;
large increments use the exact tangent. No allocation occurs in processing.
Zero depth selects the original loop without trajectory work. Boundaries reduce
depth using the maximum smoothstep derivative and knot rate, keeping the
instantaneous carrier inside the allowed positive-frequency interval. This is
not a brick-wall antialiasing guarantee: modulation has spectral sidebands.

Slow frequency wander is different: it bounds instantaneous frequency, not
accumulated phase. Existing blur uses independent per-sample phase kicks and
progressively loses coherence. The new process is neither of those.

## Tests and current gong result

- Native percussion suite: 23 tests pass, including the new movement target.
- Carrier retention and declared-energy decay tested at 32/44.1/48/96 kHz.
- Packet sharing, bounded accumulated displacement and deterministic reset tested.
- At a 15 kHz carrier, maximum depth/speed, the tested wrapped complex-sideband
  power is below −60 dB at those rates. This is a specific stress test, not proof
  of artifact-free operation at every carrier, parameter and sample rate.
- WASM configured versus prepared replay matches across 128/512-sample blocks.
- Off-state gong WAV is bit-identical to the accepted pre-change six-second WAV.
- Browser controls are visible, have tooltips and round-trip all 184 parameters.

`tools/probe_gong_motion.py` holds frequencies, levels, damping, excitation and
bloom fixed. It compares no blur, reduced blur and original blur with a grid of
movement settings on seeds 1675/1982/2586. The guarded composite score improves
from 26.71 to 24.71 for original blur plus depth 1.5, speed 100, sharing 0.25.
This is modest and does not justify calling the gong calibrated.

Standard-seed fine ridge contrast in 3–7/7–14 kHz is **5.07/4.54 dB**, versus
**3.93/3.99 dB** for the rejected blur-heavy candidate and **5.78/5.10 dB** for
the reference. Additional seeds 4903/6299/7883 retain the same general ridge
advantage. Upper centroid is still only about 4.6–4.9 kHz versus reference
6.0 kHz. The inspected difference plot still shows missing 8–12 kHz energy.
Movement alone has not fixed the upper spectral balance or fully matched its
modulation. The stronger 1.5-rad search candidate was not published; the user
requested the recommended 1-rad blend in the preset instead. The stronger
candidate's tested single/four-hit sequences invoke no browser
limiting at Master 0 dB (sample peaks −4.42/−4.16 dBFS).

Artifacts: `build/gong-bounded-motion` and `build/gong-bounded-motion-blend`.
Scripts write only diagnostic/checkpoint files; audition remains in the main
workbench. Full-bank movement adds processing cost; initial offline WASM runs
suggest roughly 25% additional render time with it enabled, but these were
concurrent search runs, not a controlled real-time benchmark.

## Subsequent listening feedback: artificial upward drift

The user hears continuous upward pitch travel rather than a stable tuned body
with later sizzle. `tools/audit_gong_bloom_separation.py` now performs read-only
ablations of the published preset at its exact standard gesture. Disabling
bounded movement barely changes the migrating band envelopes; disabling both
movement and slow wander also leaves the gross pattern. Disabling diffusion
removes the upper bloom. These checks implicate energy redistribution and its
observed spectral balance, not a monotonic pitch ramp in the movement primitive.
They do not alone reproduce or quantify the user's pitch perception.

At 44.1 kHz, with the same 4096-point STFT and 20 ms sigma smoothing of band
power, the first crossing of half each band's own maximum power in 0–2 s is:

| Band | Reference | Published gong | Movement off |
|---|---:|---:|---:|
| 0.8–2.5 kHz | 0.180 s | 0.110 s | 0.105 s |
| 2.5–5 kHz | 0.484 s | 0.299 s | 0.309 s |
| 5–9 kHz | 0.559 s | 0.614 s | 0.589 s |
| 9–14 kHz | 0.688 s | 0.858 s | 0.853 s |

These are relative **rise landmarks**, not physical excitation-onset estimates.
No normalization is applied to the audio or plotted absolute power. The middle
two bands peak about 4.6–4.8 dB above reference; the highest peaks about 4.3 dB
below. Doubling diffusion strength from 6 to 12 shortens the upper rise but raises
the 5–9 kHz peak about 14 dB above reference. Therefore neither more phase blur,
bounded movement nor an isolated diffusion-speed change repairs the temporal
layering. Reference data itself still has frequency-dependent rise; the target
is not a hard-switched or spectrally instantaneous upper layer.

The next fitting objective must explicitly protect the early 300–800 Hz pitched
body, reduce excess early mid-band prominence and bring the upper-band rise
landmarks closer together. Current frequency-neighbour diffusion couples those
tasks. Try the existing broad observation and transport controls against these
measurements before concluding a different coupling topology is necessary.
No preset or engine changes accompanied this diagnostic. Artifacts and the
inspected plot are under `build/gong-bloom-separation`.

## Literature context

[Fitz and Haken, bandwidth-enhanced sinusoidal modelling](https://www.cerlsoundgroup.org/Loris/ICMC95/BandwidthOscillators.html)
discuss carrier-plus-noise-sideband synthesis. It motivates separating carrier
retention from noise bandwidth; it does not specify our packet-sharing process
or validate this gong fit. The experiment should be judged by listening,
spectral detail, modulation and energy tests together.
