# Crash: turbulence gradient and crunchy-treble investigation

The user reported excess low-frequency noise and synthetic, almost bit-reduced
treble on both single playback and live strikes. The baseline is the
[four-point decay audition](TfPercussion-crash-decay-review.md).

## Limiter check

The saved preset's standard strike was passed through the actual browser
`lookahead_limiter_processor.mjs` in 128-frame blocks. At monitor gains −12 dB
(the visible browser setting) and 0 dB, reduction remained exactly zero. Output
was bit-identical to the scaled input after accounting for the 144-frame delay
at 48 kHz. Input sample peaks were respectively 0.12063 and 0.48025. The normal
limiter regression suite also passes. This rules out limiter distortion for
that single-hit reproduction, not every possible unsaved patch or overlapping
sequence. The limiter remains enabled.

## Sample-rate and component checks

The exact workbench Wasm engine rendered three seeds at both 48 and 96 kHz.
Comparison uses ensemble Welch power over 0.12–2 s, with equal 170.7 ms windows
in seconds. The random trajectories differ with sample rate, so a waveform
null would be invalid. The 96 kHz render contains about −73.5 dBFS integrated
power above 24 kHz in this interval. Its 8–16 kHz bands are slightly louder,
not quieter, than the 48 kHz render. Output-filter response and stochastic
realizations also differ, so this is **not** proof of alias-free DSP, but it
does not establish strong foldover as the reported fault.

Ablations disable phase diffusion, transfer diffusion or neighbour exchange
individually. Disabling phase diffusion changes upper-band energy noticeably;
disabling transfer diffusion has little effect in this preset. The recurrence
uses independent per-sample signed phase increments. Their audible texture is
still a hypothesis to test, not a diagnosed clipping/quantization bug. No DSP
algorithm was changed or claimed fixed in this pass.

## Audition change: existing turbulence controls only

The previous global turbulence was 1 with zero slope. Its upper wash was
already at the global maximum, while every resolved handle had a local
multiplier of 0.08. Therefore a global slope alone cannot greatly broaden those
upper resolved ridges. The implemented profile is

$$
T(f)=\operatorname{clamp}\left(0.8+0.18\log_2(f/2500),0,1\right)\,s_i.
$$

Here $s_i$ is the existing visible per-handle turbulence multiplier. Broad
packets retain multiplier 1. Resolved handles retain 0.08 through 1.5 kHz,
then rise linearly in log frequency to 0.35 at 4.5 kHz. These values are stored
explicitly in the modal editor/JSON, not evaluated as a new hidden macro.
The global factor is about 0.38 at 500 Hz, 0.8 at 2.5 kHz and reaches 1 near
5.4 kHz. Changing turbulence redistributes packet energy and linewidth; it is
not just an EQ, and cleaner low packets expose their central resonances more.

All modal frequencies, painted levels, T60 points, excitation, contact,
cascade, radiation and output gains remain fixed. The signed plot was inspected;
the cleaner low field also loses some already-deficient low-mid energy. The
standard spectral/attack metrics do not unanimously improve. This is an
audition of the user's requested texture direction, not an accepted fit or a
claim that crunchy treble is resolved.

Private reproducibility artifacts are in `build/crash-texture-v9/`; the selected
snapshot is `upper-0.35/`. Exact snapshot replay, the full ten-second render and
repeated strikes were checked. Only the Crash preset changes; the previous
audition remains recoverable in `build/crash-decay-v8/four-points/`.

## Quick 4x rate check

The active crash voice runs at the host rate (1x). The browser limiter's 4x
true-peak detector does not oversample synthesis. The separate oversampled FM
and self-phase-delay primitives are not in this crash's active signal path.

`build/crash-alias-v10/audition.mjs` renders the saved graded-turbulence preset
unchanged at 48 and 192 kHz: one strike, and four strikes spaced 500 ms apart,
with persistent state between strikes. Each clip lasts seven seconds. Both
observation filters remain enabled. `convert.py` converts the 192 kHz renders
to 48 kHz using a delay-compensated 1025-tap Kaiser FIR (22 kHz cutoff), verified
within 0.001 dB through 20 kHz and below -110 dB from 24 kHz. No gain matching
or normalization is applied.

Open the existing workbench with `?audition=rate-check/manifest.json` for the
four fixed-render buttons. Playback uses the normal master and limiter; these
buttons do not change the live patch/reference and do not follow UI edits.
Private generated audio and the manifest are not release assets.

Before master/limiting, single-hit peaks are -7.60 dBFS (1x) and -5.75 dBFS
(4x); repeated-hit peaks are -4.64 and -4.25 dBFS. Single-hit integrated RMS
differs by 2.46 dB, versus 0.48 dB for the repeated example. Identical controls
do not imply identical energy at different rates: filter responses and
stochastic trajectories change. Therefore this is a listening diagnostic,
not an alias-only null test or an accepted replacement preset. Production DSP
and the saved preset remain unchanged by this check.
