# Initial compact-snare fit

This note records the first calibration of `drum.snare.v1`. It is a starting
point for ear fitting, not an acceptance claim. The synthesis topology is:

```text
contact + correlated FM -> shared modal membrane -> motion-driven wire rack
          |                         |                         |
          +-------------------------+-------------------------+
                                    |
                         three-source observation -> EQ -> mono
```

The wire rack has no HIT envelope. High-passed membrane motion crosses a
threshold, follows a short engagement/release envelope, and excites correlated
tilted noise plus a normalized bank of short wire modes. Static body offset is
silent. The initial body uses a short general membrane decay and one exposed
persistent-ring frequency, decay, and level; all fitted quantities are visible
in the workbench.

## Reference and event

The first target is the local corpus cell `acoustic-snare-maple / main / v082`,
with the public event mapping strength `0.65`, location `0.30`, hardness
`0.65`, stick implement, and contact spread `0.20`. Only this cell informed the
initial parameter values. The neighbouring velocity layers remain audition
references rather than fitted evidence.

The main fitted choices are:

| Control | Initial value |
| --- | ---: |
| Membrane fundamental | 185 Hz |
| General membrane T60 | 0.12 s |
| Persistent ring | 675 Hz, 1.8 s, 0.30x |
| Direct/body velocity exponents | 2.59 / 2.23 |
| Velocity compression | 3.64x |
| Wire engagement / release | 45 ms / 30 ms |
| Wire motion threshold | 0.050 |
| Wire motion high-pass | 1,000 Hz |
| Wire range | 520--9,000 Hz |
| Wire T60 / density | 25 ms / 0.90 |
| Wire observation level | 4.00x |
| Output high cut | 10 kHz |

## Causal diagnostics

The table is measured after onset alignment and one whole-file energy match.
It deliberately keeps contact, initial decay, body, and tail separate.

| Region | Reference RMS | Synth RMS | Reference centroid | Synth centroid |
| --- | ---: | ---: | ---: | ---: |
| 0--15 ms | -17.9 dBFS | -17.0 dBFS | 279 Hz | 347 Hz |
| 15--120 ms | -26.5 dBFS | -27.6 dBFS | 1,791 Hz | 1,787 Hz |
| 120--600 ms | -50.1 dBFS | -52.4 dBFS | 679 Hz | 678 Hz |
| 600--1,500 ms | -69.7 dBFS | -71.0 dBFS | 676 Hz | 675 Hz |

The table applies a whole-file energy match of `-3.54 dB`; the workbench does
this shared reference-level match automatically while preserving the visible
`-20 dB` model-level starting control. The revised wire response now reaches
the right broad spectral region without a literal delayed burst: a causal
45 ms contact follower and squared contact law suppress the first milliseconds
and let membrane motion build wire energy. The remaining level and first-15 ms
centroid differences are explicit ear-fitting targets, not an acceptance
claim.

The separate direct and body velocity curves fix an earlier calibration bug:
the generic modal safety ceiling was active on every velocity layer and made
soft, medium, and hard hits nearly identical. The snare raises that ceiling so
it is once again a safety bound, then uses exposed velocity exponents and a
bounded saturation law to fit the event response. The neighbouring soft and
hard layers are diagnostics only; this pass intentionally calibrates just the
medium standard cell.

The reproducible development commands are:

```powershell
.\dev.ps1 build-workbench
# Run with the Node executable supplied by the configured Emscripten SDK.
& $env:EMSDK_NODE tools/render_percussion_recipe.mjs `
  build/workbench-wasm/triggerfish-percussion.mjs drum.snare.v1 `
  build/snare-initial-fit.wav sampleRate=44100 seconds=1.734 `
  strength=.65 location=.3 hardness=.65 implement=1 contactSpread=.2 seed=1575
.\.venv\Scripts\python.exe tools/inspect_audio_pair.py REFERENCE.wav `
  build/snare-initial-fit.wav
```

The live workbench is the primary listening surface. The offline renderer and
pair inspector exist to make a result reproducible and to catch causal-region
regressions; their aggregate numbers must not approve a sound automatically.
