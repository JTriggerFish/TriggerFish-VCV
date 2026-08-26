[CmdletBinding()]
param(
    [ValidateSet(
        "doctor",
        "setup",
        "build",
        "panel-preview",
        "clean",
        "dist",
        "install",
        "run",
        "smoke",
        "smoke-slop4",
        "smoke-vdpo",
        "smoke-303",
        "smoke-prog-303",
        "smoke-reverb-two-sources",
        "smoke-4072",
        "smoke-electric-piano",
        "smoke-wavefold",
        "smoke-unison",
        "smoke-scene-pack4",
        "test",
        "benchmark-er",
        "python-test",
        "shell",
        "rack-dep",
        "rack-build",
        "rack-run"
    )]
    [string]$Command = "build",

    [ValidateRange(1, 64)]
    [int]$Jobs = 8
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$repoRoot = $PSScriptRoot
$devRoot = Split-Path -Parent $repoRoot

function Get-ConfiguredPath([string]$EnvironmentVariable, [string]$DefaultPath) {
    $configuredPath = [Environment]::GetEnvironmentVariable($EnvironmentVariable)
    if ([string]::IsNullOrWhiteSpace($configuredPath)) {
        return $DefaultPath
    }
    return $configuredPath
}

$msysRoot = Get-ConfiguredPath "MSYS2_ROOT" "C:\msys64"
$rackSdk = Get-ConfiguredPath "RACK_SDK_DIR" (Join-Path $devRoot "Rack-SDK")
$rackRuntime = Get-ConfiguredPath "RACK_RUNTIME_DIR" (Join-Path $devRoot "Rack2")
$rackSource = Get-ConfiguredPath "RACK_SOURCE_DIR" (Join-Path $devRoot "Rack-src")
$msysShell = Join-Path $msysRoot "msys2_shell.cmd"
$msysBash = Join-Path $msysRoot "usr\bin\bash.exe"
function Resolve-RackExecutable {
    $configured = [Environment]::GetEnvironmentVariable("RACK_EXECUTABLE")
    if (-not [string]::IsNullOrWhiteSpace($configured)) {
        return $configured
    }

    $candidates = @()
    if (-not [string]::IsNullOrWhiteSpace($env:ProgramFiles)) {
        $candidates += Join-Path $env:ProgramFiles "VCV\Rack2Pro\Rack.exe"
    }
    if (-not [string]::IsNullOrWhiteSpace($env:LOCALAPPDATA)) {
        $candidates += Join-Path $env:LOCALAPPDATA "Programs\VCV\Rack2Pro\Rack.exe"
    }
    $candidates += Join-Path $rackRuntime "Rack.exe"

    foreach ($candidate in $candidates) {
        if (Test-Path -LiteralPath $candidate) {
            return $candidate
        }
    }
    return $candidates[-1]
}

$rackExe = Resolve-RackExecutable
$rackWorkingDirectory = Split-Path -Parent $rackExe

function Assert-Path([string]$Path, [string]$Description) {
    if (-not (Test-Path -LiteralPath $Path)) {
        throw "$Description was not found at '$Path'. See DEVELOPMENT.md."
    }
}

function ConvertTo-MsysPath([string]$Path) {
    $fullPath = [System.IO.Path]::GetFullPath($Path)
    if ($fullPath -notmatch '^([A-Za-z]):\\(.*)$') {
        throw "Cannot convert '$fullPath' to an MSYS2 path."
    }

    $drive = $Matches[1].ToLowerInvariant()
    $tail = $Matches[2].Replace('\', '/')
    return "/$drive/$tail"
}

function Invoke-Mingw([string]$Script) {
    $oldMsystem = $env:MSYSTEM
    $oldChereInvoking = $env:CHERE_INVOKING
    try {
        $env:MSYSTEM = "MINGW64"
        $env:CHERE_INVOKING = "1"
        & $msysBash --login -c $Script
        if ($LASTEXITCODE -ne 0) {
            throw "MSYS2 command failed with exit code $LASTEXITCODE."
        }
    }
    finally {
        $env:MSYSTEM = $oldMsystem
        $env:CHERE_INVOKING = $oldChereInvoking
    }
}

Assert-Path $msysShell "MSYS2 MinGW64 launcher"
Assert-Path $msysBash "MSYS2 Bash"
$repoMsys = ConvertTo-MsysPath $repoRoot
$sdkMsys = ConvertTo-MsysPath $rackSdk
$sourceMsys = ConvertTo-MsysPath $rackSource

function Invoke-PluginMake([string]$Target) {
    Assert-Path (Join-Path $rackSdk "plugin.mk") "Rack SDK"
    Invoke-Mingw "cd '$repoMsys' && RACK_DIR='$sdkMsys' make -j$Jobs $Target"
}

function Start-SmokePatch(
    [string]$Filename,
    [string]$DeviceTemplateFilename = ""
) {
    $portablePatch = Join-Path $repoRoot $Filename
    Assert-Path $portablePatch "Rack smoke-test patch"
    $localFilename = [System.IO.Path]::GetFileNameWithoutExtension($Filename) + ".local.vcv"
    $smokePatch = Join-Path $repoRoot $localFilename
    $prepareScript = Join-Path $repoRoot "tools\refresh_smoke_patch.py"
    Assert-Path $prepareScript "Smoke-patch preparation helper"
    $prepareArguments = @(
        $prepareScript,
        $portablePatch,
        $smokePatch,
        "--refresh"
    )
    if (-not [string]::IsNullOrWhiteSpace($DeviceTemplateFilename)) {
        $deviceTemplate = Join-Path $repoRoot $DeviceTemplateFilename
        if (Test-Path -LiteralPath $deviceTemplate) {
            $prepareArguments += @("--device-template", $deviceTemplate)
        }
    }
    & uv run --no-project --python 3.13 python @prepareArguments
    if ($LASTEXITCODE -ne 0) {
        throw "Smoke-patch preparation failed with exit code $LASTEXITCODE."
    }
    Write-Host "Building and installing TriggerFish..."
    Invoke-PluginMake "install"
    Write-Host "Launching Rack with $localFilename..."
    Start-Process -FilePath $rackExe -ArgumentList @($smokePatch) -WorkingDirectory $rackWorkingDirectory
}

switch ($Command) {
    "setup" {
        Push-Location $repoRoot
        try {
            & uv sync --locked --group dev --python 3.13
            if ($LASTEXITCODE -ne 0) {
                throw "Repository Python environment setup failed with exit code $LASTEXITCODE."
            }
        }
        finally {
            Pop-Location
        }
    }
    "doctor" {
        Assert-Path (Join-Path $rackSdk "plugin.mk") "Rack SDK"
        Assert-Path $rackExe "Rack executable"
        Write-Host "Rack executable: $rackExe"
        Invoke-Mingw "printf 'MSYSTEM=%s\n' `$MSYSTEM && g++ --version | head -1 && clangd --version | head -1 && gdb --version | head -1 && make --version | head -1 && cmake --version | head -1 && jq --version && python --version"
        & $rackExe --version
        & uv python find 3.13
        if (Test-Path -LiteralPath (Join-Path $rackSource "Makefile")) {
            Write-Host "Rack source: $rackSource"
        }
        else {
            Write-Host "Rack source: not cloned (plugin development is ready)"
        }
    }
    "build" { Invoke-PluginMake "" }
    "panel-preview" {
        $componentLibrary = Join-Path $rackRuntime "res\ComponentLibrary"
        $panelFont = Join-Path $rackRuntime "res\fonts\DejaVuSans.ttf"
        Assert-Path $componentLibrary "Rack component library"
        Assert-Path $panelFont "Rack panel font"
        Push-Location $repoRoot
        try {
            foreach ($moduleName in @(
                "Tf303VoiceCore", "Tf303Oscillator", "Tf4072VoiceCore",
                "TfWavefoldOscillator", "TfUnisonOscillator"
            )) {
                & uv run python tools/align_panel_labels.py `
                    --rack-runtime $rackRuntime --module $moduleName
                if ($LASTEXITCODE -ne 0) {
                    throw "$moduleName label alignment failed with exit code $LASTEXITCODE."
                }
                & uv run python tools/svg_text_to_paths.py `
                    "res-src/$moduleName.svg" "res/$moduleName.svg" `
                    --font $panelFont
                if ($LASTEXITCODE -ne 0) {
                    throw "$moduleName panel generation failed with exit code $LASTEXITCODE."
                }
            }
            & uv run python tools/svg_text_to_paths.py `
                "res-src/TfScenePack4.svg" "res/TfScenePack4.svg" `
                --font $panelFont
            if ($LASTEXITCODE -ne 0) {
                throw "TfScenePack4 panel generation failed with exit code $LASTEXITCODE."
            }
            & uv run python tools/svg_text_to_paths.py `
                "res-src/TfElectricPiano.svg" "res/TfElectricPiano.svg" `
                --font $panelFont
            if ($LASTEXITCODE -ne 0) {
                throw "TfElectricPiano panel generation failed with exit code $LASTEXITCODE."
            }
            & uv run python tools/svg_text_to_paths.py `
                "res-src/TfReverb.svg" "res/TfReverb.svg" `
                --font $panelFont
            if ($LASTEXITCODE -ne 0) {
                throw "TfReverb panel generation failed with exit code $LASTEXITCODE."
            }
            & uv run python tools/align_panel_labels.py `
                --rack-runtime $rackRuntime --module TfTransport
            if ($LASTEXITCODE -ne 0) {
                throw "TfTransport label alignment failed with exit code $LASTEXITCODE."
            }
            & uv run python tools/svg_text_to_paths.py `
                "res-src/TfTransport.svg" "res/TfTransport.svg" `
                --font $panelFont
            if ($LASTEXITCODE -ne 0) {
                throw "TfTransport panel generation failed with exit code $LASTEXITCODE."
            }
            foreach ($suffix in @("", "-30", "-38")) {
                & uv run python tools/svg_text_to_paths.py `
                    "res-src/TfProgSequencer$suffix.svg" "res/TfProgSequencer$suffix.svg" `
                    --font $panelFont
                if ($LASTEXITCODE -ne 0) {
                    throw "TfProgSequencer$suffix panel generation failed with exit code $LASTEXITCODE."
                }
            }
            & uv run python tools/render_panel_preview.py `
                --rack-runtime $rackRuntime --all `
                --documentation-directory (Join-Path $repoRoot "doc")
            if ($LASTEXITCODE -ne 0) {
                throw "Panel preview generation failed with exit code $LASTEXITCODE."
            }
        }
        finally {
            Pop-Location
        }
    }
    "clean" { Invoke-PluginMake "clean" }
    "dist" { Invoke-PluginMake "dist" }
    "install" { Invoke-PluginMake "install" }
    "run" {
        Write-Host "Building and installing TriggerFish..."
        Invoke-PluginMake "install"
        Write-Host "Launching Rack..."
        Start-Process -FilePath $rackExe -WorkingDirectory $rackWorkingDirectory
    }
    "smoke" { Start-SmokePatch "test-vdpo.vcv" }
    "smoke-slop4" { Start-SmokePatch "test-slop4.vcv" }
    "smoke-vdpo" { Start-SmokePatch "test-vdpo.vcv" }
    "smoke-303" { Start-SmokePatch "test-303-voice.vcv" }
    "smoke-prog-303" { Start-SmokePatch "test-prog-sequencer-303.vcv" }
    "smoke-reverb-two-sources" {
        Start-SmokePatch "test-room-reverb-two-sources.vcv"
    }
    "smoke-4072" { Start-SmokePatch "test-4072-voice.vcv" }
    "smoke-electric-piano" { Start-SmokePatch "test-electric-piano.vcv" }
    "smoke-wavefold" { Start-SmokePatch "test-wavefold-oscillator.vcv" }
    "smoke-unison" {
        Start-SmokePatch "test-unison-oscillator.vcv" `
            "test-wavefold-oscillator.local.vcv"
    }
    "smoke-scene-pack4" { Start-SmokePatch "test-scene-pack4.vcv" }
    "test" {
        Invoke-Mingw "cd '$repoMsys' && cmake -S . -B build/dsp-tests -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DTRIGGERFISH_BUILD_PYTHON=OFF && cmake --build build/dsp-tests -j$Jobs && ctest --test-dir build/dsp-tests --output-on-failure"
    }
    "benchmark-er" {
        Invoke-Mingw "cd '$repoMsys' && cmake -S . -B build/dsp-tests -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DTRIGGERFISH_BUILD_PYTHON=OFF && cmake --build build/dsp-tests --target triggerfish_early_reflections_benchmark -j$Jobs && ./build/dsp-tests/triggerfish_early_reflections_benchmark.exe"
    }
    "python-test" {
        $previousPath = $env:Path
        $previousGenerator = [Environment]::GetEnvironmentVariable("CMAKE_GENERATOR", "Process")
        $previousCCompiler = [Environment]::GetEnvironmentVariable("CMAKE_C_COMPILER", "Process")
        $previousCxxCompiler = [Environment]::GetEnvironmentVariable("CMAKE_CXX_COMPILER", "Process")
        $mingwBin = Join-Path $msysRoot "mingw64\bin"
        Push-Location $repoRoot
        try {
            $env:Path = "$mingwBin;$previousPath"
            $env:CMAKE_GENERATOR = "Ninja"
            $env:CMAKE_C_COMPILER = Join-Path $mingwBin "gcc.exe"
            $env:CMAKE_CXX_COMPILER = Join-Path $mingwBin "g++.exe"
            & uv sync --group dev --python 3.13 --reinstall-package triggerfish-vcv-dsp
            if ($LASTEXITCODE -ne 0) { throw "uv sync failed with exit code $LASTEXITCODE." }
            & uv run pytest
            if ($LASTEXITCODE -ne 0) { throw "pytest failed with exit code $LASTEXITCODE." }
        }
        finally {
            Pop-Location
            $env:Path = $previousPath
            [Environment]::SetEnvironmentVariable("CMAKE_GENERATOR", $previousGenerator, "Process")
            [Environment]::SetEnvironmentVariable("CMAKE_C_COMPILER", $previousCCompiler, "Process")
            [Environment]::SetEnvironmentVariable("CMAKE_CXX_COMPILER", $previousCxxCompiler, "Process")
        }
    }
    "shell" {
        $oldRackDir = $env:RACK_DIR
        try {
            $env:RACK_DIR = $sdkMsys
            Push-Location $repoRoot
            & $msysShell -defterm -here -mingw64
        }
        finally {
            Pop-Location
            $env:RACK_DIR = $oldRackDir
        }
    }
    "rack-dep" {
        Assert-Path (Join-Path $rackSource "Makefile") "Rack source"
        Invoke-Mingw "cd '$sourceMsys' && make -j$Jobs dep"
    }
    "rack-build" {
        Assert-Path (Join-Path $rackSource "Makefile") "Rack source"
        Invoke-Mingw "cd '$sourceMsys' && make -j$Jobs"
    }
    "rack-run" {
        Assert-Path (Join-Path $rackSource "Makefile") "Rack source"
        Invoke-Mingw "cd '$sourceMsys' && make run"
    }
}
