# Quiet initial highs followed by bloom

Read-only audit of saved **gong test**, snapshot
`04af4b2a-10a6-4a72-a98d-0f18f7406056`, against Dresden03. Unsaved browser
edits are not inspected. The audit below preceded the preset refinement at the
end of this note; DSP is unchanged by that refinement.

The fit has a short high-frequency onset transient, then a much deeper trough
than the reference before bloom. The engine can excite high modes initially:
the contact drives all active packets through a normalized input projection,
while diffusion subsequently redistributes stored energy into the same modes.

The audited snapshot's excitation slope is −47.94 dB/octave with a 1.5 kHz knee. Its
unnormalized gain is

$$
g(f) = (1 + (f/f_c)^2)^{s/(40\log_{10}2)}.
$$

That attenuates 5 kHz by about 86 dB and 10 kHz by about 132 dB relative to
the low-frequency limit, before contact spectrum and observation. High packets
are almost entirely dependent on later diffusion.

Changing **only initial excitation tilt to −28 dB/octave** demonstrates the
requested quiet early layer followed by bloom:

| Band | Reference, 20–100 ms | Current, 20–100 ms | Tilt −28, 20–100 ms | Tilt −28, 300–1000 ms |
|---|---:|---:|---:|---:|
| 5–9 kHz | −66.2 dBFS | −92.8 dBFS | −63.3 dBFS | −30.9 dBFS |
| 9–14 kHz | −65.4 dBFS | −113.3 dBFS | −66.7 dBFS | −36.0 dBFS |

High-band peaks remain near 0.9 seconds. This is a capability test, **not a
finished calibration**: bloom remains several dB too strong, and 2.5–5 kHz
remains weak. Normalization slightly changes the low-band level too; controls
are not fully orthogonal. Try the existing tilt before adding an engine control.

## Method

`tools/audit_gong_early_highs.py` uses the exact saved gesture and fixed gains,
rendered through WASM. Diagnostics use 10 ms disjoint RMS windows after causal
fourth-order Butterworth bandpass filters, avoiding centred-STFT look-ahead.
Filter startup transients still affect the first few milliseconds. Both signals
have a brief onset; this is not complete absence of high-frequency energy at
sample zero. The deficit is particularly clear over 20–100 ms.

Diffusion-off comparisons confirm that relaxing tilt produces direct high-mode
excitation. Enabling diffusion then gives later growth. This measures radiated
band energy, not the physical reference's internal force or energy distribution.

JSON measurements and consistently coloured Plotly curves are under
`build/gong-early-highs/`. Run with the project Python environment,
`PYTHONPATH=python`, and `EMSDK_NODE` pointing at the SDK Node executable.
No modal bars, gains, T60, presets or DSP are changed, and no listening page is
published by this diagnostic.

## Published refinement

The main workbench's **Gong — quiet onset and bloom** now uses excitation tilt
−28 dB/oct and the existing upper T60 endpoint at 1.5 seconds (previously 1.75).
All other controls, modal frequencies/levels, event, routing and gains are
unchanged from the user's `gong test`. Output EQ stays disabled. The original
snapshot remains available in the texture-trial picker.

`tools/refine_gong_early_highs.py` searches a 3×3 grid: tilt −30/−28/−26,
upper T60 1.25/1.5/1.75. The objective is RMS dB error across causal-band energy
in 20–100 ms, 100–300 ms, 300–1000 ms, 1–2 s and 2–3 s, equally weighted, for
80–800 Hz, 5–9 kHz and 9–14 kHz. It averages errors over two seeds, without gain
matching. The upper-mid scoop is measured but excluded: correcting the user's
painted response by distorting unrelated controls is not the task here.

| Band | Reference onset | Revised onset | Reference bloom | Revised bloom |
|---|---:|---:|---:|---:|
| 80–800 Hz | −17.9 | −17.6 | −25.4 | −26.0 |
| 5–9 kHz | −66.2 | −63.3 | −34.1 | −33.0 |
| 9–14 kHz | −65.4 | −66.6 | −41.4 | −39.1 |

Values are dBFS; onset is 20–100 ms, bloom is 300–1000 ms. Revised values are
the two-seed mean. Search error falls from 14.74 to 1.82 dB **on these selected
regions**, not a general perceptual score. Full eight-second renders and a third,
unoptimized seed were checked in the causal-band plots. The early layer and
subsequent growth survive the seed change, with the low body largely retained.

Limitations: 2.5–5 kHz remains too weak (about 13 dB in the bloom region), and
the revised high tail drops below the reference after about 3 seconds, where
these band levels are already around −60 to −70 dBFS. This is a targeted onset
and bloom improvement, not a claim of a finished perceptual match. It has not
been accepted by ear. Artifacts, full source snapshot, grid, held-out measurements
and Plotly curves are in `build/gong-early-refinement/`. To reproduce the original
comparison after publication, use the original `gong test` from the trial picker
as the source (`--fit exported-original.json`), not the revised default.

The actual browser limiter was also checked with one hit and four half-second
restrikes at 48 kHz: no reduction at master −12 or 0 dB. The revised preset was
loaded and exported in an isolated browser tab to verify all 184 parameters,
reference identity/gain and the exact saved gesture. The existing server serves
the rebuilt main workbench; no additional listening server was started.
