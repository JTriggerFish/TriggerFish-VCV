# Simplifying modal motion controls

The audit and reduction proposal below preceded the implemented layout trial
described at the end. The target is intuitive constructive instrument design,
not retention of every experimental parameter combination.

## What exists

Excluding packet texture, the metallic voice has ten knobs:

| Mechanism | Controls | What changes |
|---|---|---|
| Beating | Depth, rate, rate tilt | Fixed oscillator pairs and their relative excitation amplitudes |
| Pitch wander | Amount in Hz, speed | Smooth random instantaneous frequency, with accumulating phase |
| Ridge movement | Amount in radians, speed, sharing | Bounded smooth phase displacement, with a coherent ridge and sidebands at moderate depth |
| Phase blur | Amount, frequency tilt | Independent per-sample random phase increments, progressively losing coherence |

These are rotations/frequency placement, not additional envelope or energy sinks.
They can still change interference, observed loudness and nonlinear transfer.
Packet texture separately sets the surrounding oscillator population and width.

Implementation: `smooth_modal_drift.hpp`, `bounded_modal_motion.hpp`,
`modal_packet_distribution.hpp`, `crash_cymbal_parameters.cpp`, and
`stochastic_modal_field_parameters.hpp` under `src/tfdsp/percussion/`.
The workbench disables the separate legacy random-neighbour exchange.

## Deferred reduction proposal: nine parameters

This reduction is **not implemented or approved for removal**. The current trial
retains all ten controls; the user finds them useful in their clearer groups.
The proposal and its measurements remain here as a record, not a task to apply.

The user explicitly wants to retain slow detuning and is open to reorganising
the surface. That supersedes the initial seven-knob proposal: overlap is not
redundancy when slow detuning and fast shimmer must coexist.

| Group | Visible controls |
|---|---|
| Beating | Depth; speed; high-frequency speed scaling |
| Movement: Detune row | Amount in Hz; speed |
| Movement: Shimmer row | Amount; speed |
| Blur | Amount; bass/treble balance |

Retain packet texture separately. Put Detune and Shimmer together in one compact
Movement card, with aligned amount/speed rows, but keep their generators
independent. Rename the current ridge-movement experiment to Shimmer. Do not put
controls in an Advanced section or combine these distinct effects into a
many-parameter macro. Avoid an unlabelled two-dimensional pad: detune depth is
measured in Hz whereas shimmer depth is a phase excursion, so a shared vertical
axis would imply equivalence that does not exist.

Keep the separate pitch-wander path. The slow end of Shimmer covers some useful
irregular beating, but cannot substitute for independent slow frequency drift
under simultaneous fast shimmer. Nor does bounded phase displacement reproduce
the accumulated phase history of frequency wander. The actual reduction is
therefore ten parameters to nine, plus substantially less confusing layout.

Remove packet sharing from this voice: use independent modal movement directly,
not a hidden editable coefficient fixed behind the UI. Shared movement preserves
a packet's clockwork beating; independent motion disturbs those phase relations,
which is generally the behaviour wanted here. Pure regular beating remains
available with Motion off. The existing sharing mix also reduces movement
variance at intermediate settings, making amount and sharing confusingly coupled.

For the current Gong, migrate amount from 1.5 to approximately 1.295 by

$$d_{\mathrm{new}}=d_{\mathrm{old}}\sqrt{(1-c)^2+c^2},\quad c=0.15.$$

This matches marginal phase variance, not the complete joint process or sound.
It is an explicit conversion into the remaining visible amount, not a permanent
compensation coefficient. The Crash's existing detune remains unchanged; the
tested slow-shimmer substitution is evidence of partial overlap only.

Keep both spectral controls for now. Beat-rate tilt governs relative gaps between
pairs, while blur tilt governs stochastic linewidth. The Crash currently needs
opposite directions: beat tilt +0.25 and blur tilt -1.5. Tying them together would
remove an independently useful balance. A fixed beat law would reduce the count
to eight, but is less conservative for bells/gongs; do not silently substitute it.

An optional presentation improvement is a small frequency-versus-beat-speed line
with two draggable points, alongside the depth slider. This replaces the
technical rate/tilt interface without removing either degree of freedom. The
internal power law and current parameter limits should remain exact, with
numeric readouts and keyboard editing. It is an editor for two existing
parameters, not a new curve or another layer of hidden values. Initially the
compact labelled sliders are the lower-risk change.

## Tests

`tools/audit_modal_control_overlap.py` uses the actual compiled workbench WASM
with exact saved routing, controls and gestures. It writes only diagnostics to
`build/modal-control-overlap/`. Run with `PYTHONPATH=python`, the repository
Python environment, and `EMSDK_NODE`; `--only` selects one instrument or isolated
tests. No gain normalization or preset publishing.

Four-second Gong, Crash and Ride ablations use each saved seed and seed 1982.
Broad Welch bands and the existing auditory-band modulation features are
reported separately, with natural seed variation for context. These are
diagnostics, not perceptual acceptance thresholds. In particular, the existing
modulation metric starts at 2 Hz and cannot certify slow breathing equivalence.

- Gong: removing sharing and matching variance gives texture distance
  0.043–0.045, versus 0.233 between baseline seeds. Broad bands move by at most
  0.56 dB. Removing ridge movement instead gives distance 0.329–0.348.
- Crash: replacing its small wander by slow motion amount 0.6 gives distance
  0.157–0.167, versus 0.198 between baseline seeds. Broad bands move at most
  0.12 dB. This supports a listening trial, not an indistinguishability claim.
- Removing blur gives substantially larger texture changes on Crash
  (0.361–0.484) and Ride (1.14–1.20). Ride's beating controls are inactive in its
  Scattered layout, so their ablation is exactly silent; that does not make
  beating redundant for sparse paired sounds.

Additional eight-second isolated 1-kHz ring/pair tests disable energy transfer,
remove other handles and satellites, and use a known flat T60. Analysis removes
that known decay only for measuring modulation. Between seconds 1 and 7:

- Gentle bounded movement retains 93% of spectral power within ±1 Hz of the
  carrier; the tested broad blur retains 20%. The chosen strengths are not
  perceptually matched, but demonstrate separate controllable behaviours.
- A 2-Hz pair has about 99% of its measured 0.5–10 Hz envelope-modulation power
  concentrated near 2 Hz. Shared fast movement retains that regularity;
  independent fast movement lowers it to 60% while retaining pulsation depth.
- Small slow wander and slow bounded movement reduce that concentration to
  approximately 80% and 83%, respectively, with similar envelope variation.
  This demonstrates overlap, not a reason to remove independently useful slow
  detuning. The separate slow-wander controls are retained in the revised proposal.

No claim is made of matching references better, full equivalence across all
parameters, or auditory acceptance. Production DSP and saved sounds are unchanged.

## Implemented layout trial

The approved reorganisation is now in the workbench, without the proposed
parameter removal. All ten existing values remain visible and round-trip exactly;
packet sharing is labelled **Shimmer moves together** pending listening approval
of its removal. Pitch drift and Shimmer retain independent amount/speed controls.

- The left-hand control panel retains **two side-by-side control columns**:
  contact, source/output mix, T60 and final EQ on the left; size/tuning, bloom
  (excitation and energy transfer), Resonance texture, Beating, Movement
  (Pitch drift + Shimmer) and Blur on the right.
- Routing sits above both control columns, collapsed in an accordion by default.
- Kick resonance and membrane modal controls also use the right control column.
  Their excitation, tension controllers and output processing stay left.
- The separate right-hand analysis panel holds audition, the playing surface,
  and the wide modal editor and generators—not the resonance control sliders.
- Control and analysis panels scroll independently. Dependent speed/balance sliders grey out at zero
  amount without clearing their stored values. No advanced/hidden parameters.

`modal_control_presentation_tests.mjs` covers activity rules;
`modal_layout_browser_tests.mjs` covers IDs, grouping, column placement, dependent
states, the routing accordion, independent scrolling and widths 1440/1920/2560,
with screenshots.
Calibration round-trips, modal pointer interactions, the kick browser probe and
the optional WASM suite also pass. No preset JSON or production DSP was changed.

## Slider and decay-editor precision

Shimmer Amount uses squared slider travel: $d=d_{\max}p^2$ for position
$p\in[0,1]$. Its inverse is used when loading and resetting values. This gives
small amounts more room without changing the maximum, DSP units or saved fits.

The T60 graph and selected-knot slider use a softened log display coordinate
$u=\log(1+T/(1\,\mathrm{s}))$, normalized between the existing limits. About
80% of the vertical range is now available to 1–30 seconds. Short values remain
reachable at the bottom; Shift-drag still provides tenfold finer screen motion.
This is a fixed display scale, not an automatically changing zoom or a DSP
parameter. Saved seconds and the ERB/log-T60 interpolation are unchanged.

The curve is sampled from that actual interpolation before display, not drawn
as straight segments on the new nonlinear scale. The diamond still multiplies
every T60 by the same factor and stops when any knot reaches a limit; it does
not distort the decay ratios to force a parallel-looking screen shift.
