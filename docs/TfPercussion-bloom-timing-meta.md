# Bloom timing: an explicit design gesture

The workbench's **Bloom · excitation & diffusion** section groups the initial
excitation shelf with the two diffusion controls. Its timing slider is a
design-time meta-control, like Size meta, not a new DSP parameter or delay line.
Modal prominence remains in the modal editor; noisiness and damping remain
independent controls.

## Using it

The middle captures the loaded/current patch. Move left for an earlier, more
immediate body; right for a darker initial excitation and later development.
The underlying sliders update immediately. Double-click restores that captured
patch's three values. **Set centre** captures the current result without changing
the sound. Editing one of the underlying three controls also establishes a new
centre. Loading a fit captures the loaded values, without modifying them.

This is a relative design gesture, not a percentage change in measured delay.
Control-limit notices show where its motion saturates. With diffusion disabled,
the timing slider is disabled rather than silently enabling another process.

## Exact expansion

For position $p\in[-1,1]$ and the captured strength $k$, excitation shelf slope
$s$ (dB/octave) and shelf centre $c$ (Hz):

$$
k'=k\,2^{-p},\qquad s'=s-6p,\qquad c'=c\,2^{-p/4}.
$$

Each result is clamped to its ordinary descriptor range. These are the only
three changed values: `bloom_rate`, `body_brightness`,
`body_excitation_centre`. The expansion is calculated from the captured centre,
not accumulated slider deltas, so reversals and reset do not drift after clipping.

It does **not** change observation gains/bars, T60, nonlinearity, body excitation
gain, velocity, contact presentation, noisiness, or modal frequencies. There is
no loudness compensation. Saved fits contain the expanded ordinary parameters;
no timing state is required to reproduce the audio in C++/Wasm.

## Tested range and limitations

An initial four-octave strength gesture was rejected: at its late end high
energy could die before developing, making the remaining contact dominate.
The narrower one-octave gesture was checked with exact renderer sweeps using
`tools/audit_bloom_timing.py`, separately from fitting.

On the starting presets, the 3–8-kHz energy median within the first two seconds
moves monotonically across five positions: approximately 96–405 ms for crash
and 789–1022 ms for gong. These are diagnostics of those patches, not guarantees
for arbitrary instruments. A contact-dominated attack can remain at time zero
even when the body develops later. Different strikes, damping, diffusion
exponents, packet distributions and prominence can change the audible result.

Later development naturally competes with damping and can be quieter. The
gesture does not hide that trade-off by extending T60 or boosting the output.
There is no guarantee of an independent physical delay with constant bloom level.

Tests cover pure expansion, bounds, disabled transport, reversibility, visible
slider updates (including logarithmic diffusion positioning), direct-edit
rebasing, preset loading and duplicate-control avoidance. Browser checks are
silent and use disposable tabs; Wasm/browser tooling remains optional for Rack.
