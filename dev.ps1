[CmdletBinding()]
param(
    [ValidateSet(
        "doctor",
        "build",
        "clean",
        "dist",
        "install",
        "run",
        "smoke",
        "smoke-slop4",
        "smoke-vdpo",
        "test",
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
$rackExe = Join-Path $rackRuntime "Rack.exe"

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

function Start-SmokePatch([string]$Filename) {
    $portablePatch = Join-Path $repoRoot $Filename
    Assert-Path $portablePatch "Rack smoke-test patch"
    $localFilename = [System.IO.Path]::GetFileNameWithoutExtension($Filename) + ".local.vcv"
    $smokePatch = Join-Path $repoRoot $localFilename
    $refreshScript = Join-Path $repoRoot "tools\refresh_smoke_patch.py"
    Assert-Path $refreshScript "Smoke-patch refresh helper"
    & uv run --no-project --python 3.13 python $refreshScript $portablePatch $smokePatch
    if ($LASTEXITCODE -ne 0) {
        throw "Smoke-patch refresh failed with exit code $LASTEXITCODE."
    }
    Write-Host "Building and installing TriggerFish..."
    Invoke-PluginMake "install"
    Write-Host "Launching Rack with $localFilename..."
    Start-Process -FilePath $rackExe -ArgumentList @($smokePatch) -WorkingDirectory $rackRuntime
}

switch ($Command) {
    "doctor" {
        Assert-Path (Join-Path $rackSdk "plugin.mk") "Rack SDK"
        Assert-Path $rackExe "Rack runtime"
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
    "clean" { Invoke-PluginMake "clean" }
    "dist" { Invoke-PluginMake "dist" }
    "install" { Invoke-PluginMake "install" }
    "run" {
        Write-Host "Building and installing TriggerFish..."
        Invoke-PluginMake "install"
        Write-Host "Launching Rack..."
        Start-Process -FilePath $rackExe -WorkingDirectory $rackRuntime
    }
    "smoke" { Start-SmokePatch "test-vdpo.vcv" }
    "smoke-slop4" { Start-SmokePatch "test-slop4.vcv" }
    "smoke-vdpo" { Start-SmokePatch "test-vdpo.vcv" }
    "test" {
        Invoke-Mingw "cd '$repoMsys' && cmake -S . -B build/dsp-tests -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DTRIGGERFISH_BUILD_PYTHON=OFF && cmake --build build/dsp-tests -j$Jobs && ctest --test-dir build/dsp-tests --output-on-failure"
    }
    "python-test" {
        Push-Location $repoRoot
        try {
            & uv sync --group dev --python 3.13 --reinstall-package triggerfish-vcv-dsp
            if ($LASTEXITCODE -ne 0) { throw "uv sync failed with exit code $LASTEXITCODE." }
            & uv run pytest
            if ($LASTEXITCODE -ne 0) { throw "pytest failed with exit code $LASTEXITCODE." }
        }
        finally {
            Pop-Location
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
