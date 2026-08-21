from pathlib import Path
import tomllib

ROOT = Path(__file__).resolve().parents[2]


def test_windows_uv_builds_are_pinned_to_mingw_ninja():
    configuration = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    overrides = configuration["tool"]["scikit-build"]["overrides"]
    windows = [
        override
        for override in overrides
        if override.get("if", {}).get("platform-system") == "^win32"
    ]
    assert len(windows) == 1
    assert windows[0]["env"]["CMAKE_GENERATOR"] == {
        "default": "Ninja",
        "force": True,
    }
    assert windows[0]["cmake"]["toolchain-file"] == ("cmake/triggerfish-mingw.cmake")

    toolchain = (ROOT / windows[0]["cmake"]["toolchain-file"]).read_text(
        encoding="utf-8"
    )
    assert "MSYS2_ROOT" in toolchain
    assert "Visual Studio" not in toolchain
    assert 'set(TRIGGERFISH_GCC "${TRIGGERFISH_MINGW_BIN}/gcc.exe")' in toolchain
    assert 'set(TRIGGERFISH_GXX "${TRIGGERFISH_MINGW_BIN}/g++.exe")' in toolchain
    assert "CMAKE_C_COMPILER" in toolchain
    assert "CMAKE_CXX_COMPILER" in toolchain
    assert "CMAKE_CXX_COMPILER_LAUNCHER" in toolchain
    assert "CMAKE_CXX_LINKER_LAUNCHER" in toolchain
