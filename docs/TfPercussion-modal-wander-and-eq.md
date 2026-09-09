# Irregular ringing, simpler noisiness, and graphical output EQ

Current refinement: 9 September 2026. This supersedes the preceding
[slower-beating note](TfPercussion-visible-controls-and-slow-beating.md) for the
active control surface and gong/crash preset names.

## Why slower fixed pairs were insufficient

A fixed frequency difference creates a periodic beat. Lowering the difference
slows that periodicity but does not make it irregular. Increasing phase blur
instead broadens the tones into a noisy wash. The reference and candidate STFTs
show that this is not the same perceptual change.

The existing smooth-drift primitive already supplied independent, twice
continuously differentiable random frequency trajectories. But the instrument
used a percentage of frequency multiplied by the noisiness/exchange amount.
That made clean low tones barely move while allowing much wider treble movement.

The active instrument now uses **Pitch wander (Hz)** and **Wander speed
(changes/s)**. These replace, rather than supplement, the old percentage controls:

$$f_i(t)=f_i+D_i u_i(t),\qquad -1\leq u_i(t)\leq1.$$

All modes receive the same requested peak deviation $D$ in Hz. $D_i$ is only
reduced where necessary to stay between DC and the prepared frequency ceiling.
It is independent of packet noisiness, allocation, painted prominence and
phase-blur bandwidth. Each trajectory has independent targets and starting
position. Quintic interpolation joins targets with zero first and second
derivative at the joins. Speed specifies the target rate, not an LFO frequency.

Pole rotations use the existing norm-preserving Cayley construction. Wander
does not change damping radii or inject stored energy. It is free-running across
repeated strikes; only resetting/repreparing the voice resets the trajectory.
Zero wander retains the static-mode path. Generic fractional drift remains
available in `SmoothModalDrift::Prepare`; the instrument uses `PrepareHz`.

This is a **constructive approximation**, not proof that the reference contains
independent random frequency noise. Time-varying modal amplitudes, frequencies,
close partials and nonlinear coupling can all contribute. Relevant primary work:

- [Skare & Abel, DAFx 2019](https://www.dafx.de/paper-archive/2019/DAFx2019_paper_48.pdf)
  discusses crash-cymbal modal synthesis, beating and nonlinear approximations.
- [Fitz & Haken, bandwidth-enhanced sinusoidal modelling](https://www.cerlsoundgroup.org/Loris/ICMC95/BandwidthOscillators.html)
  distinguishes sinusoidal trajectories and bandwidth/noise modelling. It does
  not prescribe this particular gong wander law.

## Noisiness surface and imports

The active profile is now exactly

$$I(f)=G_{1k}(f/1000)^s L.$$

No centre parameter remains in the instrument struct, WASM descriptors, UI or
new presets. Old relaxed profiles convert with $G_{1k}=G(1000/C)^s$, preserving
the whole curve. The largest previously representable value at this pivot is
4000; the nonlinear slider gives most practical precision at small values.
Neither conversion nor DSP secretly clamps the level back to 4.

The V1 workbench surface now has 180 values after the shared-output-EQ cleanup. Import conversion is explicit and
does not mutate source snapshots. Zero legacy percentage drift converts exactly.
Nonzero percentage drift cannot be represented exactly by a single Hz depth;
such imports fail with an actionable message, not a silent approximation.
The existing blur-balance import still handles older snapshots first.

## Focused preset changes and evidence

| Preset | Peak wander | Target changes |
|---|---:|---:|
| Gong — irregular gentle ringing | 0.65 Hz | 0.5/s |
| Crash — user low ring, gentle movement | 0.3 Hz | 1/s |

Beyond the mathematically equivalent noisiness conversion, these changes touch
only wander. Painted tuning/prominence, phase blur, decay, diffusion, and gains
are unchanged. The user's crash low-frequency correction is preserved.

`tools/refine_gong_slow_beating.py --target gong --study wander-hz` screens 16
depth/speed combinations. The same experiment runs for crash. A named `--choose`
checks four seeds and creates verified candidate snapshots and exact WASM WAVs.
Artifacts are in `build/gong-wander-hz/` and `build/crash-wander-hz/`.

Across four seeds, average modulation line concentration (per-band periodicity)
changes 0.2244 → 0.1514 for gong and 0.2189 → 0.1893 for crash. Lower is less
dominated by one periodic line; this is **not** a universal perceptual quality
score. Gong reference-floor Mel changes 0.7789 → 0.7822 and decay-shape error
4.0359 → 4.1896 dB: a small regression traded for less periodic movement. Crash
Mel is essentially unchanged (1.00642 → 1.00665), as is decay error
(3.1351 → 3.1468 dB). Broad-band diagnostics can hide individual ridge behaviour.

The fixed-scale gong STFT differences were inspected. Excessively clean/strong
low-mid ridges and high-tail mismatch remain; these candidates are not certified
complete fits. Repeated quarter-note and rapid hard-strike audits remain finite.

## UI and graphical EQ

All controls remain exposed. Framed groups separate tuning/drive, packet texture,
beating/movement, and phase blur, following the existing output section style.

One shared final EQ processes the contact/body mix. It has three draggable handles: high-pass frequency,
colour frequency/gain, and low-pass frequency. Four editable numeric fields and
the enable checkbox show all existing values. Double-click resets a handle or
numeric field. No extra bands, Q, or hidden fitted parameters were introduced.
The plotted response matches the existing Butterworth cuts and Q=0.8 colour
peak. Coloured component curves and their white sum show the actual filter
shape; observation gain is separate and not included in that EQ-only curve.

The DSP is now `contact level × contact + body level × body → final EQ → model
level`. The old two per-source filter chains are removed from this voice, not
bypassed behind hidden controls. The generic observation primitives remain
available for other instruments. There is no extra delay, normalization or gain
compensation. Disabling the EQ returns the unfiltered mix.

The five `output_*` parameters replace ten per-path EQ parameters in C++, WASM,
Python and JSON. Shipped presets and imported old snapshots retain the old body
EQ as the shared curve and discard the contact EQ. This is deliberately not an
exact sound-preserving conversion; all synthesis settings and mix gains remain
unchanged. The fit measurements above predate this topology cleanup.

Regression tests compare the actual output to filtering the raw contact/body
mix at 44.1, 48 and 96 kHz, including bypass and reset. Import tests reject mixed
or incomplete old/new EQ surfaces rather than retaining dead settings.

The blue background is a **live full-mix analyser** while audio is active,
after instrument output EQ and before browser master/limiter. Both live strikes
and WAV/reference playback pass through this transparent tap. The label changes
to LIVE; when silent/suspended it falls back to the first second of the rendered
synth. Gold is the reference's first-second average, not a live reference feed.
Neither background pretends to be a separately measured contact/body input.

The analyser uses 4096 samples and the UI updates at 20 Hz. Its window affects
display response, not audio latency. Window gain and equivalent noise bandwidth
are corrected to one-sided power density; the display ceiling comes only from
the reference. No synth-dependent normalization or audio gain adjustment is
performed. Live ticks update existing histogram rectangles, not the EQ controls
or the full spectrogram. See the
[Web Audio analyser specification](https://webaudio.github.io/web-audio-api/#fft-windowing-and-smoothing-over-time)
for the Blackman window and transform convention.

## Verification and useful next UI steps

Native tests cover bounded Hz excursions, clean-mode movement, equal bass/treble
units, replay, stored energy and T60. WASM/native agreement, curve-preserving
imports, EQ response invariants and spectrum units are tested. Disposable browser
tests drag an actual EQ handle and verify values update. A muted live oscillator
test verifies the FFT follows a 1-to-4-kHz change and clears on suspension.

Suggested next improvements, not extra synth parameters:

- Replace the noisiness level/slope sliders with a small draggable frequency
  profile, retaining both numeric values.
- Add observation-only packet solo for comparing ring versus surrounding
  texture **without rebuilding allocation or interrupting energy transport**.
- Show changed-from-preset indicators and a per-section A/B/reset, so users can
  evaluate one change without losing the rest of a fit.
