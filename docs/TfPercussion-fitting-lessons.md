# Reusable fitting safeguards

These rules apply to constructive percussion models, not just the kick. They
are implemented in `python/triggerfish_percussion/` and used by the optional
developer workbench. Normal Rack builds do not depend on Python, Node or Wasm.

## First prove what the controls do

Test isolated sources, gain linearity, repeated strikes, known pitch/damping
recovery, and parameter boundaries before optimizing a real recording. Passing
these tests establishes an implementation contract, not an acoustic match.

The kick review exposed a concrete failure: switching a nearly silent modal
bar off removed a full unit of excitation normalization and raised another
mode by 3 dB. Its replacement splits continuous prominence between drive and
observation. Both weights approach zero at the displayed off floor; the modal
bank still normalizes excitation energy. The corresponding native regression
test compares -71.99 dB with off, including waveform error, not only finiteness.

## Remove redundant fitting directions

Normalized modal weights are relative. Fitting every weight permits a nearly
inaudible common scale change. Keep the loudest starting bar fixed during a
local fitting pass; fit the remaining relative levels and use the explicit
resonance gain for overall level. Do not fit master gain simultaneously with
all source gains. The UI values remain explicit and editable: fixing a value
in an optimization is not a hidden DSP coefficient.

Start damping with a shared low-dimensional curve. For the kick this is T60 at
100 Hz plus frequency slope; no per-mode multipliers. Neither output EQ nor
performance inputs are fitting variables in the current kick experiment.

## Placement before local refinement

`modal_fit_initialization.spectral_mode_candidates` proposes frequencies from
the reference's smoothed spectrum, with logarithmic-band coverage. It does
**not** claim to separate true modes from noise or measure their excitation
gains. Its long window resolves bass but smears attack timing. Use it for
starting frequencies only, then test proposals through the complete renderer.

`reference_modal_starts` generates explicit mode lists at several counts and
neutral prominence falloffs. It does not copy spectral power directly into
modal amplitudes: the source spectrum, damping and observation all affect the
measured result. Couplings in these starts are explicit unit values. The kick
tries up to 16 modes; the helper accepts a caller-selected capacity/prefix.

`Search.screen_candidates` compares discrete starts with the same fixed
reference, loss, duration and noise seeds as subsequent refinement. The current
patch is always a candidate. Candidate vectors and scores are retained. A
better initial score is not proof of the best basin: when local refinement
stalls, separately refine plausible alternatives before discarding them.
`workbench_multistart.refine_candidate_starts` implements that second step:
shortlist by raw score, then independently refit and log each alternative before
selection. `extended_modal_starts` retains the strongest current handles while
adding reference proposals, so testing upper ringing need not discard a useful
low-body fit. The kick enables this pass with `TF_KICK_MODAL_RESTARTS=1` and a
current resume file.

The previous kick start had every mode below 283 Hz, with local bounds ending
at 424 Hz. That search could never place sustained upper body ringing. Always
inspect proposed frequencies and active bounds before launching a long fit.

## Search and evidence

`Search.stage` uses bounded least squares over the exact C++/Wasm renderer.
A finite-difference sensitivity probe freezes ineffective directions; central
differences supply the Jacobian. Each stage logs effective bounds, active and
fixed values, step size, noise seeds in the saved search, objective definition,
solver status and renderer hash. Compare acceptance to the actual previous
patch, not a clamped substitute. Report the measured discrepancy by attack,
early decay and tail as well as the aggregate.

Keep absolute signal levels, reference-derived floors and colour scales fixed.
Inspect spectrograms, frequency slices, envelopes and reference/synth audio.
Do not call numerical improvement a listening-approved fit. Synthetic recovery
is necessary but insufficient: real recordings may expose missing controls,
bad starting basins or weaknesses in the loss itself.

`power_envelope.smoothed_power` evaluates the same reflected, truncated Gaussian
as the original direct convolution using FFT convolution. A measured 44.1 kHz
reference evaluation fell from 72 ms to 8 ms with maximum absolute difference
1.5e-16. Regression tests cover boundary impulses, silence and decaying signals;
this is a computational change, not a change to the loss or amplitude scales.

## Publication must fail closed

`fit_provenance.verify_candidate` requires the current recipe, Wasm hash,
descriptor set, reference identity/gain/onset and performance event. It rerenders
the saved parameter vector and independently reloads its UI fit JSON. Both must
reproduce the candidate WAV; the report reference must match the aligned source.
`workbench_fit_report.write_report` requires this verification, including in
audit-only mode. Old files merely existing is never sufficient.
`fit_publication` then pins the audio and downloadable fit to content-addressed
copies and replaces the HTML atomically. A later search checkpoint cannot change
the sound behind an already-open report's plots. Keep one writer per fitting
directory; reference audio and all these artifacts remain local.

After a model change, an old vector can be an explicitly identified warm start,
but old audio cannot be presented as current output. Never replace the user's
live preset or claim acceptance merely because an optimizer score decreased.

## Tests and entry points

- `dev.ps1 test-percussion`: native DSP, including quiet-mode removal.
- `dev.ps1 test-workbench-wasm`: native/Wasm API and JS contracts.
- `dev.ps1 test-fitting-tools`: losses, search, proposal and provenance tests.
- `TF_KICK_SELF_TEST=1` with `dev.ps1 fit-kick-start`: known C++ pitch/damping
  and coloured-noise level/decay recovery.
- `dev.ps1 fit-kick-start`: one real reference through placement and refinement.

Instrument-specific settings and measured results belong in the
[kick procedure](TfPercussion-kick-fitting.md), not in these general rules.
