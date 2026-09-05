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
| Contact | Level, width, colour, noise amount, noise T60 | Direct beater articulation and the force that excites resonance |
| Thump | Level, pitch, drop, fall time, T60 | Stable low-frequency weight with an independently timed pitch fall |
| Resonance | Level, editable modes, T60 at 100 Hz, damping slope | Persistent acoustic ringing, from absent to prominent |
| Strike/tension (advanced) | Energy pitch lift, recovery | Temporary membrane detuning proportional to stored strike energy |
| Output | EQ mode and its controls, model level | Common observation; no automatic gain matching |

Thump pitch is its settled frequency, independent of resonance frequencies.
Each of the 16 resonance slots explicitly stores frequency, relative prominence
level, and signed centre/edge strike couplings. All four are editable in the
modal panel and serialized one-to-one through the Wasm API into C++.
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
Generated spatial couplings are constructive starting weights, not measured
membrane eigenfunctions; they remain visible and editable. The same generator
is available in the metallic-body editor. It does not replace the damping curve.

Shared damping is frequency-based, independent of slot order or active count:

$$
T_{60}(f)=\operatorname{clamp}\left(
T_{100}\left(\frac{f}{100\,\mathrm{Hz}}\right)^{-s},
0.002\,\mathrm{s},30\,\mathrm{s}\right).
$$

Here $s$ is the displayed damping slope: +1 halves T60 each frequency octave;
zero is flat. There are **no individual mode decay controls or fitted
multipliers**. Output EQ is bypassed by default and excluded from kick fitting.
Optional manual radiation/multiband observation remains available.

The thump uses the existing `CorrelatedFmBurst` with zero deviation and zero
pitch jitter. Its name describes its function, not its reusable oscillator
implementation. No second FM carrier or hidden roughness control is fitted.

Contact uses the existing pulse/chirp/noise/micro-contact ingredients and
implement response. Noise T60 describes the base envelope before implement
and contact-spread shaping; the UI tooltip states that distinction.
Finite noise and thump envelopes close after an 80 dB fade, so their internal
fade duration is 4/3 of the displayed T60.

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
