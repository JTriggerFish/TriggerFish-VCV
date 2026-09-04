# Rejected compact-snare numerical start

This note records a failed numerical start for `drum.snare.v1`; it is retained
as a methodology warning, not as a calibration or current preset. The synthesis
topology is:

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

The rejected search chose:

| Control | Initial value |
| --- | ---: |
| Membrane fundamental | 185 Hz |
| General membrane T60 | 0.12 s |
| Persistent ring | 675 Hz, 1.8 s, 0.30x |
| Direct/body velocity exponents | 2.59 / 2.23 |
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

The table applies a whole-file energy match of `-3.54 dB`. Despite apparently
close RMS and centroid values, listening exposed an incoherent result: a 45 ms
wire follower delayed the defining snare response, extreme wire gain was then
hidden by a low global output, and the few aggregate descriptors did not test
the wire texture or its causal attachment to the membrane. This state was
therefore removed from the workbench defaults.

The current unverified starting point removes the wire-path delay entirely.
The bottom-head motion drives the wire contact follower directly; its 1.5 ms
engagement is the only onset smoothing. The wider 250--14,000 Hz wire response,
80 ms contact release, and 350 ms wire-mode T60 retain energy after the prompt
broadband component. The isolated 675 Hz ring remains a subordinate, separately
controllable colour rather than a delayed event.
The former saturating velocity stage has also been removed; the current
exponents may expand but cannot compress dynamics. Velocity response must be
fitted later from all three layers rather than inferred from one medium hit.

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
