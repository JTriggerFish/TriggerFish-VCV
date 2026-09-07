# Kick: contact, thump and resonance

Implemented recipe: `drum.kick.v1`. One kick replaces the separate compact/FM
and acoustic kick choices. The old CompactKick DSP primitive remains reusable,
but is disconnected from the workbench. This is a constructive instrument,
not an assertion that every acoustic kick has this exact physical topology.

## Signal flow

```text
Hit ──┬── Contact ──┬───────────────── Contact level ──┐
      │             └── Membrane ──── Resonance level ┤
      │                    ▲                           ├── Output EQ ── Master
      ├── Strike energy ───┘ (tension)                 │
      └── Thump ────────────────────── Thump level ────┘
```

The C++ voice reuses `MembraneDrum` with `DefaultKickVoiceParameters`.
It contains one contact generator per hit, one clean swept-sine thump per hit,
and a shared persistent bank of 16 membrane resonators. Eight event voices
allow overlapping contact/thump tails. There is no added runtime graph engine.

Contact's normalized body-force output drives resonance at unity. Its direct
level controls only its audible observation. Thump does not drive the membrane.
Three optional graph routes independently enable the audible contributions;
all enabled routes have unity edge gains. The three visible level controls
are their only individual volume factors. EQ and master then affect the sum.

## Controls and meaning

| Section | Controls | Meaning |
| --- | --- | --- |
| Contact | Observation choice, body drive choice, level, width, colour, noise amount, noise T60 | Independently select full/noise-only observation and full/pulse-only body drive |
| Thump | Level, pitch, drop, fall time, hold, T60, decay shape | Low-frequency weight with independent pitch and amplitude trajectories |
| Resonance | Level, editable modes, T60 at 100 Hz, damping slope | Persistent acoustic ringing, from absent to prominent |
| Strike/tension (advanced) | Energy pitch lift, recovery | Temporary membrane detuning proportional to stored strike energy |
| Output | EQ mode and its controls, model level | Common observation; no automatic gain matching |

Thump pitch is its settled frequency, independent of resonance frequencies.
Each of the 16 resonance slots explicitly stores frequency and relative prominence
level. Both are editable in the modal panel and serialized one-to-one through
the Wasm API into C++. Kick has a fixed beater: no location input or centre/edge
coefficients. The common trigger ABI ignores location for this recipe.
A level of −72 dB disables that slot's excitation and observation.
Prominence is split equally between excitation and observation: each receives
the square root of the linear bar weight minus the -72 dB off-floor. Both
therefore approach zero continuously when a mode is disabled, rather than
removing a full excitation slot at the final slider step. The bank normalizes
drive and observation separately. Bars control relative prominence; the separate
resonance level controls overall amplitude.

The modal editor's **Generate editable modes** menu supplies circular-membrane
root-ratio and harmonic-series starting layouts. Fundamental, count, level and
dB/octave falloff generate ordinary values. No formula remains active afterward.
The same generator is available in the metallic-body editor. It does not
replace the damping curve or create hidden spatial coefficients for kick.

Shared damping is frequency-based, independent of slot order or active count:

$$
T_{60}(f)=\operatorname{clamp}\left(
T_{100}\left(\frac{f}{100\,\mathrm{Hz}}\right)^{-s},
0.002\,\mathrm{s},30\,\mathrm{s}\right).
$$

Here $s$ is the displayed damping slope: +1 halves T60 each frequency octave;
zero is flat. There are **no individual mode decay controls or fitted
multipliers**. Output EQ is bypassed by default. Optional radiation/multiband
observation remains available; current diagnostic trials also test the existing
low-pass as a noise-bandwidth constraint, not a multiband correction curve.

The thump uses the existing `CorrelatedFmBurst` with zero deviation and zero
pitch jitter. Its name describes its function, not its reusable oscillator
implementation. No second FM carrier or hidden roughness control is fitted.

Contact uses the existing pulse/chirp/noise/micro-contact ingredients and
implement response. Noise T60 describes the base envelope before implement
and contact-spread shaping; the UI tooltip states that distinction.
The direct observation selector excludes pulse/chirp/grains in Noise only mode;
it does not alter their body excitation. The separate Body drive selector can
feed only the finite contact pulse to the membrane, excluding the long noise
tail from that input. It introduces no extra gain: both choices use the same
existing impulse normalization. Other drum recipes retain their full-contact
drive by default. Finite noise closes after an 80 dB fade,
so its internal fade duration is 4/3 of base T60.

After its 0.4 ms rise and optional hold, the thump amplitude is

$$
A(t)=\exp\left[-\ln(1000)\left((1-q)u+qu^2\right)\right],
\qquad u=t/T_{60},\quad 0\le q\le1.
$$

Here $t$ starts after hold and $q$ is the displayed decay shape. Zero gives an
exponential; one gives a rounded shoulder and an increasingly steep finish.
The -60 dB point stays at T60. The finite trajectory reaches -80 dB and then
closes to zero over 1 ms. This shapes amplitude, not waveform saturation.
Geometric curvature uses a recursively updated multiplier; no per-sample exp
is needed. Uncurved trajectories retain the existing float-multiply path.

## Energy, retriggers and levels

Velocity scales source amplitude linearly; there is no velocity compression.
Squared strike strength updates the separate tension state. Resonance stores
energy across hits and loses it through its modal damping. Changing resonance
prominence does not change source strength, thump, or damping.

Live parameter edits currently prepare/reset the voice, as other membrane
recipes do. **Retriggers without edits do not reset it.** This distinction
matters when testing build-up. The browser's 3 ms safety limiter is downstream
of the voice and is absent from native/plugin synthesis.

## Why these three roles

A previous same-reference ablation refitted each reduced topology instead of
just muting it. Aggregate engineering error was 4.39 dB with all three,
6.11 dB without resonance, and 4.99 dB without thump; three held-out noise
seeds supported the same ordering. Removing the old thump-to-body feed changed
the full signal by only −60.6 dB RMS. These measurements motivated this simpler
routing, not a claim of listening approval or global optimality.

See [kick fitting](TfPercussion-kick-fitting.md) for the current fitting procedure.
`dev.ps1 test-kick-architecture` repeats the ablation on the current fitted
kick. Generated reference audio and experimental artifacts remain local.
