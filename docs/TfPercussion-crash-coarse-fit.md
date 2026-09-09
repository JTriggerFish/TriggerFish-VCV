# Crash: coarse stretched-series fitting experiment

Status: **not promoted**. The saved workbench crash remains unchanged. The bloom
timing UI is implemented separately; see [its design and tests](TfPercussion-bloom-timing-meta.md).
These trials did not establish either a better crash fit or that the model
cannot fit it. No listening approval is claimed.

## Target and fixed conditions

One reference only: `crash-standard`, medium edge, private corpus A, layer v072,
take 1. Reference SHA-256:
`0446d07f06eea71428183a1e8a3507dfee8552734926e5b91c0bd4b9a10ed8c6`.
Fit the first six seconds at 48 kHz, with the existing onset alignment and
family reference gain (+42 dB). Strength is 72/127, location 1, stick implement,
hardness 0.6 and contact spread 0.2. Model level, body excitation and observation
gain stay fixed. There is no per-render normalization or velocity compression.

All nonlinear trials render the actual C++ voice through Wasm. The renderer hash
for these trials is `c1e61d2cadfc791fe796b45be1061d2d525905ec17bf1c61706408e0b83a1391`.
Checkpoints retain parameters, history, reference, event and renderer provenance.
`verify_candidate` checks current-engine reproduction before an artifact is used.

## What was fitted, and how

1. Generate whole harmonic/stretched series with the same generator as the UI:
   roots 137/175/205 Hz, 24 or 32 handles, stretch 0.25/0.5/0.75, four protected
   lower harmonics. Reject layouts that cannot fit below 15 kHz. Try both a
   neutral dynamics start and the current preset's shared controls.
2. Fit prominence through positive, piecewise-linear amplitude curves in log
   frequency, not independent bars. Compare four coordinates, then six at
   120/400/1500/4000/8000/15000 Hz. All local noisiness multipliers remain one.
3. Fit shared diffusion strength/exponent, excitation slope/centre, packet
   noisiness/slope/spread/bandwidth and two T60 endpoints using bounded
   least-squares with measured finite differences. These controls really render
   through the engine; there is no neural surrogate voice.
4. Refine the low-root family (125.24/128 Hz, 24 handles, stretch 0.8–0.95).
   The retained simple trial is **128 Hz, 24 handles, stretch 0.9, core 4**.
   Its frequencies stay exactly on that generated series throughout refinement.
5. Test contact tone/noise balance, width, tilt and direct level, then sparse
   shared damping: first endpoints only, then one 2500-Hz interior knot, then
   two at 500/2500 Hz. Compare bounded L-BFGS-B and Powell on the exact scalar
   objective. No per-mode decay fitting.
6. Refine the broad prominence curve using Auraloss Mel multi-resolution STFT
   with a 1-dB band/time bloom guard. Autograd is used only on the exact affine
   observation basis, verified against actual renders. This is a spectral loss,
   **not a learned perceptual embedding or a listening test**.

Search ranking is Mel loss + 0.05 times the norm of `SpectralBloomLoss`, whose
band/time definition is saved with each run. Ranking alone cannot approve a fit:
separately inspect attack, band decay, reference/model/difference plots and
repeated hits. Training includes seed 1396978464; final prominence refinement
also uses +101. Audit +307 and +911 are additional, unused fitting seeds.

## Results and rejection

For the standard seed, the old preset's full Mel loss is 1.3805. The preferred
simple three-point-T60 trial is 1.4771. Its 0–400-ms Mel loss worsens from 1.0835
to 1.3966, and independent band-decay shape error from 3.7748 to 4.2530 dB.
A composite balance diagnostic improves from 10.9608 to 10.4530, illustrating
why accepting a fit from one aggregate improvement is unsafe.

The simple trial's difference and decay plots were inspected. They still show
missing early/mid development around 1.5–6 kHz and an overly persistent 6–16-kHz
tail, alongside overly prominent low harmonics. Repeated hits remain finite,
but the raw rapid-hard sequence peaks at +6.61 dBFS before the browser's output
protection. These are not release-quality calibration results.

Further rejected experiments:

- Three small measured low/mid centre corrections: worse than the exact series;
  discarded rather than escalating to independent mode fitting.
- Moving a broad prominence knot: no useful improvement.
- Joint six-amplitude/sparse-T60 fitting: only a small improvement in the
  band/time objective; did not resolve the main envelope errors.
- Lower global packet noisiness and narrower spreads: worse matching. The
  hypothesis that a simple texture reduction would fix the missing midrange
  was not supported by this screen.

The two-interior-T60 trial improves the combined two-seed score only about
0.0017 over one interior point, insufficient reason to prefer the extra knot.
None of these results justifies replacing the workbench preset.

## Reproduction and artifacts

Use the optional development environment and current Wasm build via `dev.ps1`;
set `EMSDK_NODE` to its Node executable and keep BLAS/Torch at one thread.
These tools never publish presets, start servers or change Rack builds.

```powershell
python tools/fit_coarse_metal.py crash --output build/crash-coarse-warm
python tools/fit_coarse_metal.py crash --source build/crash-coarse-warm/candidate --roots 125.24 128 --counts 24 --stretches .8 .85 .9 .95 --output build/crash-coarse-low-root
python tools/refine_coarse_metal.py crash build/crash-coarse-low-root/candidate --output build/crash-coarse-finish
```

The completed Powell continuation is in `build/crash-coarse-powell/decay-1`;
its parent contact checkpoint is `build/crash-coarse-finish/contact`. Logs and
`search.json` record the actual staged runs, including the resumed bounds fix.
`difference.png`, `decay.png` and `refit-audit.json` are in that decay checkpoint.
The baseline is `build/crash-coarse-warm/baseline`.

`refine_coarse_geometry.py`, `polish_coarse_curve.py`,
`fit_joint_coarse_decay.py` and `screen_coarse_texture.py` are diagnostic
experiments, not extra mandatory fitting stages. Their `candidate` means best
within that experiment, **not better than the incoming or published preset**.
Keep the baseline and review before promotion. Joint damping variants rebuild
interior knots from the endpoints; they do not preserve an existing interior curve.

The next useful investigation is a controlled band-development test of the
current excitation/transport/damping path, separating contact and body output.
The present experiments do not identify whether the remaining error is mainly
search conditioning, the restricted family, or transport behaviour. Adding
many individually fitted modes would obscure that distinction.
