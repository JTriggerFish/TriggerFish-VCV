# Workbench edit feedback

The measurements below describe the UI/progressive-analysis pass. A subsequent
[modal-motion optimisation](TfPercussion-modal-motion-optimization.md) roughly
halved the current gong's DSP time, to about 1.25 seconds for the complete render.

Offline synthesis already runs in a dedicated Web Worker, independently of the
live AudioWorklet. It runs as fast as the DSP permits, not in real-time blocks.

Measured in a disposable Brave tab with the current 7.7-second Gong target:

| Stage | Before | Current progressive rendering |
|---|---:|---:|
| Control settling delay | 220 ms | 60 ms |
| First analysed image after a single edit | about 2.7 s | about 110 ms |
| Complete offline DSP render | about 2.4 s | about 2.5 s |
| STFT | about 50 ms for full signal | typically a few ms per new segment |
| Canvas drawing | about 10 ms | about 5–10 ms |

These are local spot measurements, not portable performance guarantees. The
change improves time to visual feedback, **not underlying DSP throughput**.
Continuous changes are coalesced until a 60 ms pause; obsolete synthesis stops
at the next worker yield. No queue of outdated full tails is rendered.

`offline_render.mjs` runs the same C++ DSP in chunks, first publishing roughly
0.125 seconds of audio (at least one complete selected FFT window) and
subsequently a growing prefix at about 120 ms wall-time
intervals. The worker yields through MessageChannel to service newer edits,
without the repeated-timer minimum-delay penalty. Previews are display-only;
the completed audio is used for the saved rendered snapshot. Chunked output is
tested against an uninterrupted WASM render and is bit-identical.

`ProgressiveStft` caches completed frames by render identity and FFT settings.
Only newly available frames are calculated. Incomplete right-edge windows are
omitted until more audio arrives. The final frame values and peak are tested
against an uninterrupted STFT for all three window types.

The display retains the preceding spectrogram behind the new prefix, including
when a new edit interrupts a partial render. A grey write edge separates current
data from the older tail in every applicable viewing mode. This display history
is never used as analysis input or saved audio. A complete render replaces it
entirely, including when the new requested duration is shorter. Unrendered time
without any history stays blank rather than being coloured as a deficit. The
reference colour ceiling remains fixed.

The difference palette uses black for exact equality, amber for reference
excess, cyan for synthesis excess, with a fixed ±24 dB saturation range and a
labelled sign legend. Its 257-entry palette has an exact black midpoint.

Reproduce with `tools/profile_workbench_updates.mjs [optional-canvas.png]` using
the existing workbench server and disposable debugging browser on port 9223.
It tests one edit and a burst of twenty edits without touching the playing tab
or generating audible sound. Run the test suite through
`./dev.ps1 test-workbench-wasm`.

## Remaining DSP cost

Review correction: snapshots taken while a render is pending save their current
controls but do not attach PCM from the previous completed render. Restoring a
cached snapshot cancels any older in-flight render, and Save fit exports the
visible edits even after selecting an earlier snapshot. The browser calibration
probe covers these cases. Display history remains separate from saved audio.

`tools/profile_gong_render_components.mjs` performs diagnostic ablations without
changing a saved fit. Median of three warmed two-second WASM renders at 44.1 kHz:

| Diagnostic configuration | Render time |
|---|---:|
| Current gong | 655 ms |
| Ridge movement bypassed | 223 ms |
| Diffusion bypassed | 535 ms |
| Both bypassed | 99 ms |

These are not alternative sound settings or acceptable shortcuts. They identify
bounded per-mode movement as the largest additional cost in this case; they are
not a measured speedup available while preserving the sound. No C++ DSP was
changed in this UI pass. Full seven-second-plus renders still take approximately
2.5–2.9 seconds locally. The subsequent optimisation linked above preserves the
tested audio while substantially reducing the movement/propagation cost.
