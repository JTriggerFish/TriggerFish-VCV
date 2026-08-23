# Development

TriggerFish Elements targets VCV Rack 2 and currently builds against the Rack
2.6.6 SDK. Plugin builds must use the compiler and ABI expected by the SDK for
their operating system.

The standalone DSP tests use CMake, Ninja, and C++17. Python 3.13 and the
optional analysis bindings are managed by [uv](https://docs.astral.sh/uv/).

## Windows

Rack plugins for Windows must be compiled with MSYS2 MinGW64, not MSVC. The
PowerShell helper expects the following sibling directories by default:

```text
development-directory/
|-- Rack-SDK/
|-- Rack2/          # Optional Rack runtime for run and smoke
|-- Rack-src/       # Optional Rack source checkout
`-- TriggerFish-VCV/
```

MSYS2 uses its standard installation directory by default. These locations can
be overridden with `MSYS2_ROOT`, `RACK_SDK_DIR`, `RACK_RUNTIME_DIR`, and
`RACK_SOURCE_DIR`. The `run` and `smoke-*` commands prefer an installed Rack
Pro executable, fall back to `Rack2/Rack.exe`, and can be directed to a
specific edition with `RACK_EXECUTABLE`. Keep the repository and SDK paths free
of spaces because the Rack Makefile toolchain does not support spaces in
absolute paths.

Run the usual tasks from PowerShell:

```powershell
.\dev.ps1 doctor       # Verify tools and configured paths
.\dev.ps1 build        # Build plugin.dll
.\dev.ps1 panel-preview # Regenerate and render editable module panels
.\dev.ps1 install      # Install into the Rack user plugin directory
.\dev.ps1 run          # Install and launch Rack
.\dev.ps1 smoke-slop4  # Install and open the Slop4/four-VCO patch
.\dev.ps1 smoke-vdpo   # Install and open the two-VDPO patch
.\dev.ps1 smoke-303    # Install and open the sequenced 303 voice patch
.\dev.ps1 smoke-prog-303 # Replace Foundry with Prog Sequencer in the 303 patch
.\dev.ps1 smoke-reverb # Install and open the 303-through-room-reverb patch
.\dev.ps1 smoke-4072   # Install and open the MIDI-playable 4072 voice patch
.\dev.ps1 smoke-wavefold # Install and open the Wavefold Oscillator patch
.\dev.ps1 smoke-unison # Install and open the Unison Oscillator patch
.\dev.ps1 smoke-scene-pack4 # Install and open the Scene Pack 4 patch
.\dev.ps1 dist         # Build the release .vcvplugin package
.\dev.ps1 test         # Build and run standalone C++ DSP tests
.\dev.ps1 benchmark-er # Benchmark generated-FIR ER at 1, 4, and 8 sources
.\dev.ps1 python-test  # Build Python bindings and run pytest
.\dev.ps1 clean
.\dev.ps1 shell        # Open a configured MinGW64 shell
```

The default command is `build`. Use `-Jobs 4`, for example, to change build
parallelism.

VS Code's tasks call this PowerShell helper. The checked-in editor settings are
platform-neutral; ensure MinGW64 clangd is available to the editor on Windows.
Running `dev.ps1 test` once generates `build/dsp-tests/compile_commands.json`
for clangd.

The optional Rack source commands are `rack-dep`, `rack-build`, and `rack-run`.
Normal plugin development only requires the SDK.

## Linux

Install a C++17 compiler, GNU Make, CMake, Ninja, curl, jq, tar, zstd, and uv.
Download the Rack 2.6.6 Linux SDK and export its absolute location:

```bash
export RACK_DIR=/path/to/Rack-SDK
```

Build, package, and install the plugin with the Rack SDK Makefile workflow:

```bash
make -j"$(nproc)"
make dist
make install
```

Run the standalone test suites with:

```bash
cmake -S . -B build/dsp-tests -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON \
  -DTRIGGERFISH_BUILD_PYTHON=OFF
cmake --build build/dsp-tests --parallel
ctest --test-dir build/dsp-tests --output-on-failure

uv sync --locked --group dev --python 3.13
uv run pytest
```

Compare the current VDPO integrator with the legacy BDF implementation using:

```bash
uv run python tests/python/benchmark_vdpo.py
```

Benchmark the Slop linear-Hz pitch conversion using:

```bash
uv run python tests/python/benchmark_slop.py
```

Benchmark the 2x/4x diode-ladder implementations and their high-drive quality
difference using:

```bash
uv run python tests/python/benchmark_diode_ladder.py
```

Compare the reduced BA662 VCA with the pre-existing TriggerFish transistor and
OTA VCA cores, with output level matched before calculating the residual:

```bash
uv run python tests/python/benchmark_ba662.py
uv run python tests/python/benchmark_tb303_oscillator.py
```

Run the slower, manual transistor-level BA662-clone comparison with ngspice:

```bash
uv run python tests/python/reference_ba662_spice.py --ngspice /path/to/ngspice
```

This downloads the official Nexperia matched-transistor models only into the
ignored `build/ba662-reference` directory. It is not part of the plugin build
or required for the automated tests.

The GitHub Actions workflow performs both a Linux Rack SDK package build and
the standalone tests, so Linux compatibility is checked continuously.

## Smoke-test patches

[test-slop4.vcv](test-slop4.vcv) is a MIDI-controlled, enveloped Slop4 voice
using four Fundamental VCOs. Two VCOs are at unison and two are one octave
higher.

[test-vdpo.vcv](test-vdpo.vcv) contains two MIDI-controlled, enveloped
Slop/VDPO/VCA voices. One VDPO is self-resonating. The other is forced by a
Fundamental VCO; click the Push module to cycle its forcing waveform through
sine, saw, and square.

[test-303-voice.vcv](test-303-voice.vcv) combines 303 Oscillator
(`Tf303Oscillator`) and 303 Voice Core (`Tf303VoiceCore`) with Impromptu
Modular's Clocked and Foundry. Its 16-bar song repeats the Gibber Acid motif
over i, -iv, and -V for 8, 4, and 4 bars. Separate Fundamental LFOs provide
slow cutoff and resonance cycles; a 14 Hz LFO is connected to
linear filter FM with its attenuverter at zero. Install Impromptu Modular from
the Rack Library before opening it.

[test-prog-sequencer-303.vcv](test-prog-sequencer-303.vcv) keeps the tuned
Clocked, oscillator, voice, modulation, Room Reverb, and stereo output path but
replaces Foundry with Prog Sequencer. Its visible 16-line program is both a
playable test and a compact-syntax regression fixture.

[test-4072-voice.vcv](test-4072-voice.vcv) connects a Fundamental saw
oscillator to 4072 Voice Core, with MIDI pitch tracking both oscillator and
filter and MIDI gate driving the two internal envelopes.

[test-wavefold-oscillator.vcv](test-wavefold-oscillator.vcv) is a MIDI-playable
Wavefold Oscillator patch exposing the oscillator and folded outputs together.

[test-unison-oscillator.vcv](test-unison-oscillator.vcv) is a MIDI-playable
Unison Oscillator patch exposing its mono and stereo mixes.

[test-scene-pack4.vcv](test-scene-pack4.vcv) feeds four fixed-pitch Fundamental
oscillators into Scene Pack 4 and exposes its packed AUDIO and X outputs on a
scope while monitoring the first packed audio channel in stereo.

The Slop4, VDPO, 4072, Wavefold, Unison, and Scene Pack patches use Rack Core,
Fundamental, and TriggerFish modules; the 303 patch also uses Impromptu Modular.
Their final stereo masters are set to -6 dB. `dev.ps1 smoke` remains an alias
for `dev.ps1 smoke-vdpo`.

The checked-in patches do not name a MIDI controller, audio interface, or
driver. Select the appropriate MIDI and audio devices after opening one. Keep
monitor volume low on first use.

The patch files are generated from `tools/generate_test_patches.py`. Edit that
source and run `uv run python tools/generate_test_patches.py` when changing the
patch topology or defaults.

Editable panels are in `res-src`. Rack's NanoSVG renderer does
not support SVG text, so regenerate the runtime asset after changing labels,
using the DejaVu Sans font bundled with a Rack runtime (or the system copy on
Linux):

```text
uv run python tools/svg_text_to_paths.py \
  res-src/Tf303VoiceCore.svg res/Tf303VoiceCore.svg \
  --font /path/to/DejaVuSans.ttf
```

An optional `--bold-font` can supply real bold outlines; otherwise the helper
adds a small outline to bold labels. The font itself is not copied or committed.

Render the panel together with the actual Rack knobs, switch, trimmers, ports,
and screws without launching Rack:

```text
uv run python tools/render_panel_preview.py --rack-runtime /path/to/Rack2
```

Pass `--module Tf303Oscillator` to select one module, or use `--all` to render
the complete collection. The tool reads widget positions from each module
source, embeds the installed Rack component artwork into a self-contained SVG,
and writes SVG and PNG previews under the ignored `build/panel-preview`
directory. `--documentation-directory doc` also refreshes the stable PNG names
used by the README. Set `PANEL_PREVIEW_BROWSER` if Edge, Chrome, or Chromium is
not found automatically. On Windows, `dev.ps1 panel-preview` regenerates both
editable runtime panels and all documentation previews in one command.

The smoke commands use ignored `*.local.vcv` launch copies so machine-specific
state never enters Git. Every command atomically reconstructs its local copy
from the selected checked-in fixture, ensuring it cannot deliberately launch
stale topology or program text. Only audio and MIDI device selections are
carried forward, and only when the old local patch is readable JSON; parameter
changes and extra local modules are intentionally discarded.

## Vendored sequencer parser

TfProgSequencer uses cpp-peglib from `vendor/cpp-peglib`. The exact upstream
release, commit, file hashes, and update procedure are recorded in
`vendor/cpp-peglib/README.triggerfish.md`. No system parser package is needed.
After an update, run `./dev.ps1 test` and `./dev.ps1 dist`; the distribution
must contain `vendor/cpp-peglib/LICENSE`.

Future MIDI performance binding, instrument interpretation, richer pattern
branches, CV expressions, tuning, and related realtime constraints are retained
in the internal
[Prog Sequencer development roadmap](docs/TfProgSequencer-planned-features.md).

## Python and pre-commit

Install the development environment and Git hooks once per clone:

```text
uv sync --locked --group dev --python 3.13
uv run pre-commit install
```

Black formats Python and notebook code. nbstripout removes notebook outputs,
execution counts, and transient metadata before commit. Run all checks with:

```text
uv run pre-commit run --all-files
```

For interactive analysis with SciPy and Matplotlib:

```text
uv sync --locked --group dev --extra analysis --python 3.13
uv run python
```
