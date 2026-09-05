# Percussion gain and energy contract

Current contract for metallic, membrane, snare and FM-kick recipes. This
supersedes older references to calibrated output boosts, contact-port
projections, normalized painted observation, and state-energy caps. Units are
constructive synthesis conventions, not SI mechanical units.

## Source ports

Direct radiation is audio. Body drive is modal excitation per sample. Both
contain the same sources at the same declared amplitudes; only their unit
conversion differs. Conversions are prepared once per event, never measured
from a playing tail or fitted separately for a random seed.

| Source | Body conversion |
| --- | --- |
| Half-sine | Reciprocal discrete envelope sum: amplitude is total impulse |
| Chirp | Reciprocal damped envelope sum, not signed waveform sum |
| White-noise burst | Reciprocal square root of envelope squared sum and uniform-noise variance |
| Micro-contact burst | Stochastic-window normalization, undoing the direct-audio carrier normalization |
| Membrane FM feed | Reciprocal peak-normalized amplitude-envelope area in samples |
| Continuous wire modal feed | Reciprocal square root of sample rate: noise excitation power per second |

The chirp normalization retains cancellation and resonant frequency selectivity.
Noise tilt and micro-contact occupancy still shape delivered excitation; they
are not compensated by seed-dependent AGC. Uniform bipolar noise has variance
$1/3$. Micro-contact direct audio uses analytical low-pass variance normalization.
Its brightness cutoff is 800–20,000 Hz, Nyquist constrained, rather than
implicitly changing with sample rate. Tilt is a complementary shelf with gains
$10^{-t/40}$ and $10^{t/40}$; automatic RMS makeup was removed.

Spatial excitation satisfies $\sum_m |q_m|^2=1$. This does **not** guarantee equal
energy from every contact waveform. For propagated state $\tilde z$ and added
excitation $d$, actual work is

$$E_{new}-E_{propagated}=2\operatorname{Re}\langle\tilde z,d\rangle+\|d\|^2.$$

The cross term remains: restrikes can reinforce or oppose existing motion.
Body excitation is amplitude gain: 4x gives 16x isolated linear energy with
other controls held fixed. Cascade activation uses the unit input-vector energy
as its reference. Cascade/exchange remain passive; T60 and mute remove energy.

## Observation and density

Metallic bars are literal observation amplitudes, not normalized makeup. Every
state in a packet observes with its handle's gain; the input weights already
share its excitation budget. With decorrelated phases, expected power is
proportional to $\sum_m g_m^2|q_m|^2$. Refining a packet does not multiply that
budget. Coherent peaks can differ: density is neither AGC nor peak normalization.
Scaling every bar scales output, not excitation. An exactly silent handle is
deactivated, as indicated by the editor.

Membrane and wire spatial projections use unit norms, without the previous
0.25/2 and 1.25 calibration factors. No stored-energy ceiling rejects subsequent
force. Valid high-energy states remain linear and decay passively. Non-finite
and subnormal sanitization is numerical fault handling, not musical saturation.

## Gain inventory

| Stage | Visible owner | Meaning |
| --- | --- | --- |
| Source character | Implement, hardness, impact balance/width, advanced contact controls | Deterministic recipe laws for mix, duration and spectrum |
| Incident amplitude | Velocity, body excitation; membrane contact/FM body levels | Amplitude plus intentional timbral response |
| Body transport | Cascade rate and energy acceleration | Passive redistribution, no gain recovery |
| Observation | Painted bars, body/contact/wire levels, EQ | Literal gains; EQ may boost |
| Voice output | Model level | Exactly $10^{L/20}$; no metallic x174 |
| Comparison | Editable Reference gain | Fixed corpus default, shared across velocity layers; never applied to the model |
| Listening | Master | Before the final limiter; not instrument state |
| Protection | Browser limiter | Absent from native/Rack model and analysis renders |

Removed: eight contact-port coefficients, extra implement body-drive multiplier,
brush stiffness makeup, tilt RMS makeup, metallic x174, membrane force/output
boosts, wire output boost and sensitivity divisor, reference-selection matching,
and membrane/wire energy ceilings. Wire sensitivity now has literal units (old
UI values divided by eight); the divisor is not retained in processing.

Fixed performance laws remain algorithm definitions, not independent hidden
preset state. `CrashCymbal::ContactParameters` derives source settings from
implement/hardness/velocity; membrane and kick recipes do likewise. Numerical
filter/window constants and source-family mappings have not been replaced by
dozens of independent fit knobs. Changing those laws changes the recipe and
must not be disguised as level fitting.

Native raw-parameter fitting is a diagnostic surface, not proof that a fit
round-trips through the smaller workbench macro surface. Published fits must be
verified by rendering their actual UI/JSON parameter values.

## Browser and comparisons

```text
reference file → displayed fixed corpus calibration → master → limiter → device
model output ──────────────────────────────────────→ master → limiter → device
```

Lookahead is 3 ms total, including the 32-frame reconstruction analysis delay.
Protection uses a sliding maximum, linked channel gain, 100 ms dB-domain recovery
and 30 ms hold (hold adds no audio latency). A 4x sinc detector estimates
intersample peaks. The nominal -1 dB ceiling reserves the $\cos(\pi/8)$
peak-sampling margin plus 0.1 dB interpolation tolerance. This is documented
protection headroom, not model makeup or formal BS.1770 certification. The UI
shows pre-limiter peak and gain reduction, warning above 3 dB reduction.
Severe limiting is not a valid fitting reference.

Reference calibration is applied once on import and included in both plots and
audition. It is recorded alongside the original source hash in snapshots. The
same gain is used for every velocity layer. A corpus's gain never amplifies the
model: browser testing caught 42 dB of reference compensation previously being
applied to the synth too. That shared-trim design was removed.

There is no level matching, automatic or one-shot. Changing samples, layers,
velocity or timbre never adjusts Model level. Reference gain is independently
editable and remembered across a corpus's layers; double-click restores its
fixed corpus default. Edits always rescale original samples, not the previous
scaled buffer. Snapshots store both visible gains alongside the source hash.
Spectrogram caches include reference gain so explicit gain edits update the
reference correctly. Model edits never renormalize the reference colour scale.

Previous fits are not compatible calibrations after this unit change. Metallic
and snare/tom factory model levels are -6 dB; acoustic and FM kick are -12 dB.
Metallic body excitation/observation default to unity. Gong uses visible 4x body
observation and 0 dB output; its contact observation was raised by the same
factor. Snare's three observation gains were raised together, retaining their
balance and leaving wire excitation unchanged. These levels were reset without
changing frequencies or decay and **without claiming a refit**. User snapshots
are not deleted or silently migrated.

The fixed reference-start level audit (48 kHz, the selected standard cell;
eight full-strength hits spaced 250 ms apart for the restrike column) gives:

| Start | Single peak, dBFS | Restrike peak, dBFS |
| --- | ---: | ---: |
| Crash | -13.6 | -4.1 |
| Snare | -15.4 | -8.8 |
| Acoustic kick | -10.3 | -4.9 |
| Gong | -16.5 | -8.4 |
| Ride | -19.2 | -12.0 |
| Hi-hat | -14.7 | -9.4 |

These are raw model peaks, before the -12 dB browser master or limiter, not
perceptual loudness matching or a guarantee for arbitrary patches. The plain
tom factory start peaks at -14.9 dBFS for a strength-0.8 strike. Reproduce with
`./dev.ps1 audit-workbench-starts`; the audit does not alter any parameters.

## Regression tests

- Half-sine impulse across duration and 44.1/48/96/192 kHz.
- Noise/micro-contact ensemble response at a fixed physical modal frequency.
- Single 80 Hz mode remains sinusoidal at 4x; energy scales by sixteen.
- Packet refinement preserves expected incoherent power and literal bar gains.
- Membrane linearity above its former ceiling, with passive free decay.
- Existing repeated-hit, velocity, routing and prepared-state tests.
- Browser unity transparency, impulses, sustained bass/high-band overload up to
  amplitude one million, bass purity, and independent 8x reconstructed peaks.
- Native/Wasm agreement and browser reference-selection/level-reset checks.

Subsequent fitting must use raw renders and explicit gains, never limiter output,
per-hit normalization or a hidden restoration of the old boosts.
