# Visible controls and slower gong beating

Follow-up: [modal wander and graphical EQ](TfPercussion-modal-wander-and-eq.md)
implements the noisiness reduction proposed here and supersedes the preset/UI details.

## Scope — 9 September 2026

All active instrument sound sections now start open. Nested “Advanced” sections
were removed, including contact and radiation. Modal controls are grouped into
tuning/drive, packet texture, beating, phase blur and pitch wandering. Parameters
remain one-to-one with the patch and DSP; this pass does not silently replace
editable values with constants. Inapplicable paired controls remain visible but
disabled. Users can still collapse a top-level section themselves.

The crash and ride presets are unchanged in this pass.

## Which controls really overlap?

There is a genuine redundancy in the current relaxed noisiness profile:

$$I(f)=G\left(\frac{f}{C}\right)^s L.$$

Here $G$ is Packet noisiness, $C$ is Noisiness centre, $s$ is slope, and $L$
is the local multiplier. For a fixed pivot $P$, the identical curve is

$$G'=G\left(\frac{P}{C}\right)^s,\qquad I(f)=G'\left(\frac{f}{P}\right)^sL.$$

**Recommended next parameter reduction:** replace level + centre + slope with
“Noisiness at 1 kHz” + slope. This removes a parameter, rather than hiding it.
It needs an explicit snapshot/JSON conversion and review of the new level range:
some legal old curves require a converted level above the current maximum of 4.
Do not silently clamp them. This conversion is not implemented in this UI pass.

Keep these distinct:

- Noisiness changes centre-versus-satellite energy; spread changes frequency
  support; density changes how many oscillators fill it. They interact but are
  not interchangeable. Local allocation allows different packet densities.
- Beat rate sets a pair's frequency separation; depth sets partner balance;
  rate tilt changes how separation varies across packets.
- Blur changes phase coherence; blur tilt shapes that effect across frequency.
- Smooth drift moves frequencies along continuous random trajectories. It is
  not phase blur. Both current gong and crash use zero drift; it is a candidate
  for an explicit future removal if listening tests find it unnecessary, not a
  hidden “advanced” setting.

No additional macro or knob was added.

## Gong diagnosis and bounded retune

The preceding fit reduced shared periodicity but left the base doublet rate at
3 Hz, with a +0.25 frequency tilt. That was not a test of the user's slower-rate
request. It gives 3, 4.24, 6 and 8.49 Hz at 125, 500, 2000 and 8000 Hz.

The new preset **Gong — slower gentle beating** changes exactly two sound values:

| Parameter | Previous | Current |
|---|---:|---:|
| Beat rate at 125 Hz | 3 Hz | 0.65 Hz |
| Beat depth | 0.5 | 0.3 |

The +0.25 tilt is retained: rates are now 0.65, 0.92, 1.30 and 1.84 Hz at those
frequencies. Pitches, mode levels, allocation, blur, diffusion, damping and all
gains remain unchanged. Other nearby oscillators can still beat at other rates.

## Reproducible checks and limitations

`tools/refine_gong_slow_beating.py` archives the starting parameters, then screens
24 combinations: rate 0.35/0.65/1/1.5 Hz, depth 0.15/0.3/0.5, tilt 0/0.25. It uses
six-second actual workbench WASM renders at the standard reference gesture and
fixed playback gains. `--choose rate-0.65-depth-0.3-tilt-0.25` checks four seeds
and writes verifiable before/candidate checkpoints. Nothing is auto-published.

`modulation_signature` now reports separate relative-envelope depths in 0.5–3 Hz
and 3–12 Hz. Tests use known 0.75 Hz versus 4 Hz AM tones to ensure these are not
treated as equivalent. The four-second analysis region has 0.25 Hz resolution;
sub-0.5 Hz movement and smooth bloom curvature are not identified as individual
beats. Broad-band envelopes can obscure motion in individual ridges.

Across four seeds the chosen candidate reduces mean rapid modulation in every
strike. Reference-fixed Mel error is 0.77997 → 0.77892 and decay-shape error is
4.2216 → 4.0359 dB. However **not every modulation metric improves**: absolute
fast-depth mismatch increases from 0.0502 to 0.0693, partly because the reference
has strong motion near 3 Hz in 180–320 Hz. Shared-line concentration also rises
slightly as the deliberate beats move into fewer slow bins. Neither is a reason
to call the whole result a perfect fit or blindly minimize one scalar.

The fixed-scale STFT difference was inspected (8192 samples at 44.1 kHz,
1024-sample hop). Excessively straight low/mid ridges and high-tail mismatch
remain. This is a slower, gentler audition candidate, not listening approval.

Generated artifacts live under `build/gong-slow-beating/`: baseline parameters,
screen/validation tables, candidate/reference WAVs, repeated-hit audit, and
`candidate/difference.png`. Audition uses the main workbench, not another server.
