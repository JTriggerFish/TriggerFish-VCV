# Crash and gong: fitting the relaxed modal field

This is a parameter-fitting experiment on the existing C++ instrument, not a
new sound generator. Only the crash and gong factory fits are in scope. The
user's saved ride and the other instruments stay unchanged.

## Fixed experiment

Each target uses its existing reference, onset, source gain, strike velocity,
location, implement, hardness and random seed. Crash uses 10 seconds at 48 kHz;
gong uses 8 seconds at 44.1 kHz, including zero padding after the source ends.
There is no candidate normalization, output gain fitting or limiter in the
comparison WAV. The browser's existing master and safety limiter remain active.

The complete parameter vector, reference identity, event, Wasm hash and search
history are stored at each checkpoint. `verify_candidate` checks that both the
saved JSON and saved WAV reproduce with the current engine before publication.

## What is fitted, and how

1. Enable relaxed turbulence. Screen 32 texture candidates, then alternate
   bounded finite-difference least squares for texture, excitation/transport and
   shared damping with observation-amplitude optimization.
2. Refine using the actual **auraloss mel multi-resolution STFT objective**:
   spectral convergence plus log-magnitude error at FFT sizes 512, 2048 and 8192.
   Bounded Powell searches operate on exact Wasm renders for the nonlinear
   texture, transport, damping and contact controls. This is a scalar objective,
   not a one-element residual passed to a high-dimensional least-squares solver.
3. Between those stages, optimize positive painted-bar observation amplitudes
   with PyTorch analysis gradients. The audio basis consists of actual C++
   renders; a mixed-amplitude render and directional gradient check validate it.
   This differentiates the observation mixture and analysis, **not the nonlinear
   instrument**. Every final proposal is rerendered through Wasm.
4. Rebalance across two training seeds, including shared damping and energy
   acceleration, then constrain mel observation polishing. The existing guard
   allows at most 0.5 dB additional RMS error in each band's absolute energy and
   relative decay, and each early attack bin, separately for each training seed.
   This is a non-regression constraint relative to the balanced patch, not a
   claim that the balanced patch already matches the reference. It does not
   reshape the playback envelope or normalize the audio.

Texture means global turbulence, slope, packet spread, phase bandwidth and local
exchange. Transport means upward cascade rate, transfer diffusion and initial
excitation tilt. Contact means ping/noise balance, width, noise tilt and ping
pitch. Existing active T60 knots can move vertically; no knots are added and no
per-mode decay multipliers are fitted. Painted frequencies and local turbulence
multipliers stay fixed in these passes.

The first screening/least-squares experiment varied turbulence centre as well
as level. That contains a redundant direction in the relaxed mapping:

$$
I(f)=L\left(\frac{f}{f_c}\right)^s t_{\mathrm{local}}.
$$

At fixed slope, scaling the centre by $a$ and the level by $a^s$ changes nothing.
Subsequent refinement holds the selected centre fixed and varies level/slope.
The centre remains an editable UI control; this is a fitting-coordinate choice,
not a hidden playback coefficient.

## Review, not a score-only acceptance test

The original spectral/envelope composite and the library mel objective disagree
on some candidates. In particular, the first gong composite improvement made
mel error worse; mel-only bar optimization then reduced too much high-frequency
bloom. Neither result alone establishes a good fit.

Inspect reference/candidate spectrograms on the same scale, their difference,
absolute band envelopes, and raw 0–1/1–3/3–10/10–30/30–100 ms attack energies.
Check other stochastic seeds and quarter-note/rapid repeated hits. Background
noise at the end of the reference is not necessarily cymbal decay. No automatic
rule turns an improved score into listening approval.

Useful partial improvements may replace the main workbench preset while later
checkpoints are being explored. The frozen previous patches remain in the local
search artifacts. A new seed inspected during refinement is a validation seed,
not thereafter an untouched test sample.

## Reproduce

Use `dev.ps1` to prepare/build the optional workbench and its existing server.
Normal Rack builds do not require these Python analysis tools.

```powershell
$env:EMSDK_NODE = 'D:/dev/emsdk/node/24.19.0_64bit/node.exe'
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
.venv/Scripts/python.exe tools/refit_relaxed_metals.py crash --output build/relaxed-refit-v13/crash --rounds 2
.venv/Scripts/python.exe tools/refine_metal_perceptual.py crash --resume build/relaxed-refit-v13/crash/pass-2 --output build/relaxed-refit-v14/crash --rounds 2
```

Replace `crash` with `gong` for the other target. Use a fresh output directory for
another experiment. `plot_spectral_difference.py` and `plot_decay_comparison.py`
produce the comparison figures from a complete checkpoint directory.
`review_metal_refit.py` records the fixed baseline, extra seeds and repeated hits.
These tools do not publish presets or launch additional web servers.

`polish_metal_guarded.py` performs the final two-seed balance/guarded-mel pass.
It takes a complete `--resume` checkpoint, a fresh `--output` directory, and an
explicit `--seed-offset` for the second training seed. Any previously audited
seed used here is thereafter training data and must be labelled accordingly.

See [the relaxed mapping](TfPercussion-relaxed-turbulence-experiment.md),
[the spectral/envelope objective](TfPercussion-metallic-fitting-v2.md) and
[the analysis toolkit](TfPercussion-analysis-toolkit.md) for implementation and
measurement details.
