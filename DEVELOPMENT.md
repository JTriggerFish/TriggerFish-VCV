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
`RACK_SOURCE_DIR`. Keep the repository and SDK paths free of spaces because the
Rack Makefile toolchain does not support spaces in absolute paths.

Run the usual tasks from PowerShell:

```powershell
.\dev.ps1 doctor       # Verify tools and configured paths
.\dev.ps1 build        # Build plugin.dll
.\dev.ps1 install      # Install into the Rack user plugin directory
.\dev.ps1 run          # Install and launch Rack
.\dev.ps1 smoke        # Install and open the smoke-test patch
.\dev.ps1 dist         # Build the release .vcvplugin package
.\dev.ps1 test         # Build and run standalone C++ DSP tests
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

The GitHub Actions workflow performs both a Linux Rack SDK package build and
the standalone tests, so Linux compatibility is checked continuously.

## Smoke-test patch

[test.vcv](test.vcv) contains every TriggerFish module and uses only Rack Core
and Fundamental dependencies. It provides two MIDI-controlled, enveloped
voices: Slop/VDPO/VCA and Slop4/four VCOs/VCA. Two VCOs are at unison and two
are one octave higher. A final stereo master is set to -6 dB.

The checked-in patch does not name a MIDI controller, audio interface, or
driver. Select the appropriate MIDI and audio devices after opening it. Keep
monitor volume low on first use.

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
