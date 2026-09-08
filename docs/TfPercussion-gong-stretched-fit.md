# Gong: stretched-harmonic starting fit

**Historical experiment:** this run used the earlier power-law stretch, before
the protected harmonic core was added to the generator. The saved JSON still
contains those explicit fitted frequencies. It has not been silently regenerated.
See the [current generator law](TfPercussion-spectral-diffusion.md).

This was an **audition candidate**, not an approved calibration. The main
workbench now uses the [protected-core coarse fit](TfPercussion-gong-coarse-fit.md).
The preceding fit is archived in `build/gong-stretched-grid/baseline`.
No DSP topology, output gain, velocity mapping or reference gain was changed.

## Target and initialization

The target is **Gong Dresden 03**, using the first six seconds of its mono
downmix at 44.1 kHz and the existing −6 dB source-family gain. Reference SHA-256:
`a36721cff1f77c22484ac026330aa4da953826ff0e92376d7fb756e27d945147`.
The fixed gesture is strength .76, location .55, hardness .35, mallet .5,
contact spread .3. There is no per-file normalization.

`fit_stretched_gong.py` calls the **actual workbench JS generator** through the
renderer bridge. It does not reproduce the formula independently in Python:

$$ f_n=f_0 n^{1+s}. $$

It screened roots 55, 82.4069, 110 and 123.4708 Hz; counts 16, 24 and 32;
and upper endpoints 10 and 14 kHz. Stretch follows from the count and endpoints.
Initial bar levels interpolate the previous fit in log frequency; all local
noisiness multipliers start at one. The best three starts receive an actual
observation-level fit before choosing the layout.

The selected start has **24 handles, root 123.4708 Hz, stretch .488586**.
Root/stretch refinement subsequently produces 123.9681 Hz / .48139,
ending at 13.739 kHz. A strictly stretched series misses important low ridges,
so this is an initialization, not a constraint imposed on the final instrument.

## What was fitted, and how

1. Fit shared excitation tilt/centre, diffusion strength/energy dependence,
   packet noisiness/slope/spread/phase bandwidth and the **two existing T60
   endpoints**. Bounded least-squares uses finite differences of actual Wasm
   renders, with logarithmic coordinates for wide positive ranges. No local
   damping multipliers or additional damping knots are fitted.
2. Alternate with observation-bar fitting. An affine basis of actual renders
   gives exact STFT cross-power matrices, including interference. The cached
   objective/Jacobian is checked against direct renders. These levels do not
   affect stored energy or transport.
3. Fit root/stretch jointly. Early 400-ms and full-duration library Mel losses
   alone proved insufficient: they sometimes rewarded sharpening the low note.
   Adding linear-bin MR-STFT improved discrimination but did not remove this
   trade-off. Those sharp-note candidates were **not published**.
4. Inspect the reference attack spectrum (16,384-sample Hann Welch window,
   65,536 FFT, first 400 ms). Place low centres at **121.8, 343.9 and 553.8 Hz**;
   reuse two quiet upper handles for **247 and 374.1 Hz**. The 247-Hz ridge is
   clearer later in the sample. Keep the handle count at 24 and refit levels
   and shared dynamics after this geometric change. This is explicit,
   target-specific ridge selection, not an automatic decomposition claim.
5. Balance the four strong attack ridges in ±6-Hz measurement windows. Update
   their visible observation bars by the measured dB error, at most 4 dB per
   step; retain a step only when the measured ridge error decreases. This fixes
   excess 344-Hz ringing and deficient fundamental/374-Hz energy without EQ.
6. Freeze the observation bars below 700 Hz and refine upper levels using
   auraloss Mel over seeds 1675 and 1776, with a 1-dB per-cell bloom guard.
   Frequencies, global dynamics, damping and gains stay fixed in this last pass.

All resulting values are ordinary editable UI/JSON parameters. The generator
settings are **not hidden runtime controls**; replacing modes in the generator
will replace the edited fit, not reconstruct it exactly.

## Checks and remaining mismatch

Final checkpoint: `build/gong-stretched-reviewed`. Current-renderer baseline
and candidate use the same Wasm binary, source, duration and audio levels.
Both fixed-scale spectrogram/difference plots and band-envelope plots were
inspected. Four synthetic seeds and repeated hits were rendered.

| Diagnostic (lower is better) | Previous fit | Candidate |
|---|---:|---:|
| Full Mel, standard seed | 1.1019 | 1.0917 |
| Full Mel, four-seed mean | 1.2343 | 1.2398 |
| Attack linear-bin MR-STFT, four-seed mean | 2.0218 | 1.9825 |
| Attack Mel, four-seed mean | 1.2382 | 1.2425 |
| Bloom/envelope residual norm, standard seed | 6.1965 | 6.3454 |

This is a **low-ridge improvement, not an overall numerical win**. Attack
linear-bin MR-STFT improves on all four seeds; full-spectrum results remain
mixed. The broad low decay and delayed upper bloom are reasonably aligned.
The midrange tail is still too weak/short, and the upper attack has excess
energy before the main bloom. Fine texture and individual ridge trajectories
remain different. Numerical measures do not substitute for listening approval.

The T60 endpoints are approximately 6.917 and 2.866 seconds. The global
diffusion strength is 5.158, energy dependence .223, noisiness .796 and slope
.646. These are fit results, not universal gong settings.

Single-hit peaks across the four seeds range from −11.65 to −8.33 dBFS.
Eight quarter-note hits peak at −6.55 dBFS; eight rapid full-strength hits
peak at **+2.62 dBFS before master/limiter**. The float output is finite, but
this is not a no-limiting headroom pass at unity master gain. No hidden gain
reduction was inserted to conceal it. Strength .25/.5/.75/1 checks show rising
upper-band energy; these are not additional recorded-velocity calibrations.

## Reproduction

The commands below describe the historical run. The current generator/fitter
now uses the protected-core law and will not reproduce these old starting grids.
Use the archived checkpoint parameters to rerender this candidate exactly.

Build with `dev.ps1 build-workbench`; Python/Torch remain optional development
dependencies. Point `EMSDK_NODE` at the configured Emscripten Node executable.
The local reference server is needed by the exact renderer bridge.

The checkpoint chain is: `fit_stretched_gong.py` →
`refine_stretched_gong_pitch.py` → `polish_metal_bloom_perceptual.py` →
`refine_stretched_gong_pitch.py --fine-ridges` →
`refine_stretched_gong_attack.py` → `polish_metal_bloom_perceptual.py` →
`balance_gong_attack_ridges.py` →
`polish_metal_bloom_perceptual.py --seeds 1675 1776 --fixed-below 700`.
Each stage saves full parameters, objective settings and its parent path.

Use `review_metal_refit.py`, `audit_metal_strength.py`,
`plot_decay_comparison.py` and `plot_spectral_difference.py` on the resulting
checkpoint. Browser audition remains in the main workbench, not separate pages.

The reusable addition is `AttackRidgeLoss`: published auraloss MR-STFT with
linear bins and 4096/8192/16384 windows on the first 400 ms. Its identity,
small pitch error, fixed gain, time crop and invalid-input checks are tested.
It is a diagnostic/complement, **not a reliable stand-alone pitch acceptance
test**. Observation fitting can now freeze explicitly measured bars, with tests
that inactive bars remain inactive and only valid observation keys are frozen.

## Subsequent fitting direction

The independent upper-bar refinements above supplied too many weakly identified
degrees of freedom. They are exploratory history, **not the recommended default
for the next fit**. Start from a low harmonic core with progressively stretched
upper modes. Fit individually only ridges that can actually be followed in the
reference. Treat the dense/noisy upper region as a broad spectral distribution:
coarse level/tilt or a few smooth frequency bands, shared damping, noise amount
and bloom development. Do not interpret stochastic STFT detail as measured
individual oscillator frequencies or amplitudes. A dedicated coarse upper-fit
parameterization is now implemented and exercised in the subsequent
[coarse fit](TfPercussion-gong-coarse-fit.md).
