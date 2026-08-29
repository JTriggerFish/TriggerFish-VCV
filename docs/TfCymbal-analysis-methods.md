# Cymbal analysis and fitting methods

This is a research and requirements note for the replacement analysis system.
The experimental cymbal analysis and fitting scripts have been removed; none of
their estimates are accepted targets. The methods below must be reimplemented
as small modules and validated independently before they constrain the ride.

## Decision

Use a hybrid analysis rather than one universal ridge detector:

1. Model the first contact transient separately.
2. Estimate resolvable low- and mid-frequency components as exponentially
   damped sinusoids, including close pairs that produce beating.
3. Stop interpreting the dense upper spectrum as individually measurable
   modes.  Compare its spectral envelope, density, modulation, and decay on a
   perceptual filterbank instead.
4. Validate both parts with resynthesis and held-out cymbals.  A parameter fit
   is not accepted merely because it lowers one spectrogram error.

This division is physically and perceptually appropriate.  Measurements have
found more than 100 modes plus split partners in an 18-inch cymbal, while
nonlinear cymbal models and measurements describe a wide-band energy cascade.
At high modal overlap, individual peaks in an FFT cease to be reliable objects.

## Modal estimator candidates

### First baseline: SAMPLE and BeatsDROP

[SAMPLE](https://github.com/LIMUNIMI/SAMPLE) is an MIT-licensed Python package
designed specifically to estimate frequency, amplitude, and decay from impact
sounds.  Its tracker rejects increasing trajectories, can join interrupted
tracks, and fits decay with a hinge regression so that the recording noise
floor does not flatten the estimate.  Its BeatsDROP extension can split a
single amplitude/frequency trajectory into a close modal pair when beating is
present.

It is a good rapid baseline, not ground truth.  Version 2.2.0 currently assumes
NumPy 1.x (`np.reshape(..., newshape=...)`) and needs a small compatibility
patch for this project's NumPy 2.x environment.  It should remain an offline
analysis dependency rather than part of the Rack plugin.

A local synthetic smoke test used modes at 347, 350.5, 711, and 1234 Hz with
known independent decays plus decaying broadband noise.  With a 16384-sample
analysis window, ordinary SAMPLE collapsed the close pair to 347.01 Hz;
BeatsDROP recovered 346.30 and 349.80 Hz and their two decay constants closely.
The same long window missed the deliberately short 1234 Hz mode.  This is a
useful warning: cymbal analysis must be multiresolution and sub-band, even when
using a purpose-built tracker.

### Reference estimator: sub-band ESPRIT

Sirdey et al., [Modal analysis of impact sounds with ESPRIT in Gabor
transforms](https://www.dafx.de/paper-archive/2011/Papers/61_e.pdf), model an
impact response directly as exponentially damped sinusoids plus noise.  They
apply ESPRIT in Gabor sub-bands so the long decay can be observed at manageable
cost and the noise is closer to white within each band.  Their metal-sound
example extracted 246 candidates and retained 152 after audibility pruning,
with a perceptually close resynthesis.

Ege, Boutillon, and David's [high-resolution modal
analysis](https://arxiv.org/abs/0909.0885) combines sub-band conditioning,
ESTER model-order selection, and ESPRIT.  It separates plate modes in cases
where ordinary Fourier analysis performs poorly and reports useful results up
to roughly 70 percent modal overlap.  This is the stronger reference method for
checking whether SAMPLE has merged or invented low/mid modes.

The open [Modal-estimation](https://github.com/orchidas/Modal-estimation)
toolbox contains frequency-zoomed and frequency-warped ESPRIT, but is MATLAB
code and does not state a licence.  It is useful as a research reference; a
small NumPy/SciPy implementation or a clearly licensed alternative is safer
for a reproducible project tool.

### Diagnostic ridge views

- `librosa.reassigned_spectrogram` sharpens a spectrogram by reassigning bins
  to local instantaneous time and frequency.  It is useful for visualizing
  close tracks and frequency drift, and for seeding candidates, but it is not
  by itself a modal estimator.
- [ssqueezepy](https://github.com/OverLordGoldDragon/ssqueezepy) provides
  MIT-licensed STFT/CWT synchrosqueezing and ridge extraction.  It is useful for
  detecting velocity-dependent frequency drift and nonlinear energy transfer.
  Its AM/FM-component assumption is least trustworthy in the dense cymbal wash.
- [loristrck](https://github.com/gesellkammer/loristrck) wraps the Loris
  reassigned partial tracker and exposes frequency, amplitude, phase, and
  bandwidth.  It is useful for deterministic/noise decomposition and
  resynthesis experiments, but is GPL and is not needed in the shipped plugin.

These views should corroborate the damped-mode estimate, not replace
high-resolution parameter identification.

## Recommended measurement pipeline

### Pre-processing

- Keep the original level whenever the source has not normalized each hit.
- Onset-align from a high-frequency energy derivative, then retain a short
  pre-onset noise-floor region.  Do not align on the absolute waveform peak;
  hard edge hits can bloom after the physical contact.
- Analyze close, overhead, and room microphones independently.  A mode seen in
  several channels/takes is more credible; concatenating channels creates
  false modes.
- Separate approximately 0-8 ms contact, 8-50 ms transition, and the later
  resonant tail.  Exact boundaries should be selected from the data.

### Resolvable modes

- Run at least three time/frequency resolutions and overlapping frequency
  bands.  Short windows recover fast high modes; long windows estimate close,
  slowly decaying modes.
- Use SAMPLE plus BeatsDROP first.  Use sub-band ESPRIT/ESTER on selected
  low/mid bands as the independent reference.
- Merge duplicate band/window estimates in log-frequency space.  Retain a mode
  only if it has enough reconstructed energy and is stable across adjacent
  windows, repeated strikes, or microphones.
- Store frequency, initial amplitude, phase where meaningful, decay constant,
  frequency drift, estimator uncertainty, and evidence count.  Do not reduce a
  cymbal to a list of frequencies alone.
- Determine the modal/statistical crossover from increasing modal overlap and
  reconstruction residual, rather than hard-coding one frequency.  A likely
  starting search range is 3-8 kHz, varying by cymbal, velocity, and decay time.

### Dense residual and decay

Subtract only the accepted deterministic resynthesis, then describe the
residual in an energy-conserving ERB or fractional-octave filterbank:

- energy share and spectral quantiles versus time;
- EDT and robust piecewise decay slopes in each band, including the noise-floor
  intersection and confidence interval;
- spectral flatness and crest factor;
- temporal modulation/burst rate and inter-channel coherence;
- residual autocorrelation and late-tail diffuseness;
- velocity-dependent spectral centroid, bandwidth, and upward energy transfer.

This is the part that should constrain the dispersion/dense-response model. Peak count
per kHz is too sensitive to FFT size and threshold to serve as a target.

### Comparing model and reference

Use several losses with explicit roles:

1. Match mode sets one-to-one using Hungarian assignment on log-frequency,
   then compare frequency in cents, decay in log-seconds, and energy in dB.
   Penalize unmatched modes and mode-count errors.  The official 2026 DAFx
   modal-estimation challenge uses the same basic idea: log-frequency optimal
   assignment followed by frequency, decay, gain, and mode-count errors.
2. Compare the dense residual with multiresolution log-STFT and mel/ERB losses.
   [auraloss](https://github.com/csteinmetz1/auraloss) provides Apache-licensed
   multiresolution STFT and perceptually weighted mel-STFT implementations.
3. Compare band-energy decay curves directly.  A spectrogram score can look
   good while all high bands decay too quickly, which is the audible failure we
   have already encountered.
4. Compare onset/contact features separately so a broad hard-stick transient
   cannot compensate numerically for a sparse, over-ringing tail.
5. Report medians and distributions across repeated hits, not only a single
   aggregate score.

Before any gradient fit, loudness-align only according to the experiment being
run.  Preserve absolute velocity-to-level curves when comparing dynamics; use a
small gain-only alignment when comparing timbre.  Never independently
normalize every time-frequency frame.

## Dataset use and generalization

Bitwig is part of the fitting corpus, alongside Iowa, Salamander, VCSL, and any
controlled Toontrack renders that are in scope.  It is not merely a final
listening check and it is not the sole target.

Fit shared trends on multiple independent sources and reserve whole cymbals,
not random individual hits, as validation data.  For example, train size and
location trends on three Bitwig rides plus Salamander and validate on the
fourth Bitwig ride and Iowa.  Rotate the held-out instruments.  This prevents
the model from learning one library's microphone, processing, or velocity
mapping.

The fitted representation should distinguish:

- cymbal-level structure: size/profile family, stable mode distribution, and
  frequency-dependent damping;
- strike location: modal amplitudes and contact spectrum, with frequency shifts
  permitted only where repeated measurements support them;
- velocity: contact energy and bandwidth plus measured nonlinear frequency
  drift/energy transfer;
- implement: contact duration and spectral envelope, not a global pitch shift.

## Validation before fitting the real cymbals

Build synthetic impacts with known ground truth and vary close-pair spacing,
decay, SNR, colored wash, frequency drift, onset duration, and room tail.  Score
mode precision/recall, cents error, decay error, energy error, and residual
reconstruction.  The [2026 DAFx Parameter Estimation Challenge reference
implementation](https://github.com/LOGUNIVPM/1st-DAFx-Challenge) supplies a
current metal-plate generator and an official evaluation script, so it is a
particularly useful external benchmark.

Only after an estimator passes those tests should its ridges become tuning
targets.  Real-sample resynthesis and blind listening remain the final check,
because a real cymbal is nonlinear and its upper tail is not exactly a fixed
bank of independent sinusoids.

The acceptance rule is an AND, never an average: the salience-weighted modal
gate, the broadband spectrotemporal gate, and the matched listening gate must
all pass independently. Better flatness, density, or broad-band T20 cannot
compensate for the wrong dominant ridges. Likewise, matching a sparse set of
mode frequencies cannot compensate for the wrong wash or time evolution.

## Implementation order

1. Implement synthetic estimator fixtures with known damped modes, close pairs,
   drift, coloured dense response, and noise floors.
2. Implement and validate onset/segmentation and ERB-band energy/decay descriptors.
3. Add one damped-mode baseline and an independent sub-band ESPRIT reference.
4. Require reconstruction, uncertainty, and cross-velocity evidence before an
   estimate becomes a target.
5. Add named multi-objective comparisons and deterministic real-versus-model
   auditions only after the forward components can be rendered in isolation.

## Nonlinear wash experiments

Static linear density is not sufficient. Controlled 64-line VFM candidates and
a 4,096-mode upper-field candidate could move individual decomposition metrics
in the right direction, but neither reproduced the joint flatness, phase
predictability, attack, and band-envelope statistics of the Iowa ride. In
particular, adding static FDN coupling made the 2-8 kHz poles more stable and
more exposed. Increasing the explicit mode count reproduced stationary-energy
share but did not supply the measured finite bandwidth or upward spectral bloom.

This agrees with nonlinear thin-shell literature:

- Ducceschi and Touzé model the cymbal/gong regime as nonlinear modal energy
  transfer and wave turbulence, not independent damped sinusoids.
- Cirio et al. split a small nonlinear low-frequency simulation from a cheaper
  high-frequency model and synthesize turbulent detail with a phenomenological
  upward diffusion of spectral energy.
- Skare and Abel report high-frequency components emerging roughly 100-300 ms
  after a crash impact and propose velocity-dependent modal coupling as a
  real-time approximation.
- Conan et al. give a stable coupled-filter formulation in which only modal
  energy above a threshold is redistributed through a normalized matrix. Random
  phase on positive transfer terms evokes plate wave turbulence without making
  the persistent linear modes themselves stochastic.

Any replacement nonlinear experiment must expose its energy transfer and
intermediate audio explicitly. Soft or late motion should return toward the
stable body response, while energetic early motion may transfer energy upward
and decorrelate receiving high-frequency components. This is a hypothesis to
test, not a retained implementation.
