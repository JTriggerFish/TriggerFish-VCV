# Metallic workbench review

This checkpoint collects the spectral-diffusion experiment, saved metallic
presets, fitting diagnostics, harmonic generator and error/decay-editor work.
It does not introduce another fit or a new bloom macro.

## Review fixes

- Removed tiny-weight denominator floors from packet coordinates and silent
  cell allocation. Very dark excitation must not move a packet toward DC or
  destroy energy arriving from another packet. Preparation uses double
  precision before normalizing back to float state weights.
- Quiet positive stored energy is rescaled rather than ignored or re-seeded;
  a zero target clears the state. Regression tests compare ordinary and tiny
  excitation/state amplitudes under the same linear diffusion law.
- Worker failures release pending render/analysis bookkeeping. Decode errors
  use event listeners; assigning `onmessageerror` did not receive the injected
  event in the tested browser. Tests drop a response, inject an error, and
  verify that the next edit renders successfully.
- Superseded reference loads no longer report a false failure or restore an
  obsolete fitted gesture after rapid preset changes.
- Connected the diagnostic-audition test to the optional Wasm suite, corrected
  its error-object assertion, added a finite-gradient guard to coarse Mel
  refinement, and removed a redundant descriptor wrapper.

No output gain, normalization, velocity response or saved fit parameter was
changed by these review fixes. The final gong render is finite, peaks at
−8.51 dBFS before presentation, and differs from its previous six-second render
by 0.0744% relative waveform L2. This is a numerical regression check, not
listening approval. Previous fit reports retain their original renderer hashes;
they must be regenerated before use as scores for the rebuilt DSP.

## Verification

Through `dev.ps1`: 22 percussion tests, 530 Python tests, two native workbench
tests, twelve Wasm/JavaScript tests, and native/Wasm signature parity.
Silent disposable-browser tests cover all four metallic presets, rapid preset
switching, worker recovery, generator limits, protected harmonics, noisiness,
persistent errors and fine T60 dragging. No user-tab reload or audio autoplay.
Python, browser tests and Wasm remain optional for normal Rack builds.

## Next control-surface discussion

The difficulty separating high-mode prominence from bloom timing remains.
Consider an explicit design-time bloom meta-control that writes visible
excitation and diffusion parameters, while leaving observation prominence
independent. It must not insert a delayed burst, hide runtime coefficients or
claim an exact delay independent of strike strength. This is a proposal to
test, not an implemented control.

The [current gong fit](TfPercussion-gong-coarse-fit.md) intentionally prioritizes
a simpler protected series; its reported errors are **not better** than the
previous detailed fit. The experiment still needs listening and further work.
