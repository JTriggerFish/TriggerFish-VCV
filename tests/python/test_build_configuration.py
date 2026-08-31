from pathlib import Path
import re
import subprocess
import sys
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


def test_macos_10_9_sources_avoid_unavailable_libcxx_entry_points():
    # VCV's MacOSX12.3 SDK annotates these C++17 libc++ entry points as
    # requiring macOS 10.13 even though Rack plugins deploy back to 10.9.
    unavailable = (
        ("std::any", r"\bstd::any\b"),
        ("std::any_cast", r"\bstd::any_cast\b"),
        ("std::visit", r"\bstd::visit\b"),
    )
    offenders = []
    for pattern in ("*.cpp", "*.h", "*.hpp"):
        for path in (ROOT / "src").rglob(pattern):
            source = path.read_text(encoding="utf-8")
            for name, pattern in unavailable:
                if re.search(pattern, source):
                    offenders.append(f"{path.relative_to(ROOT)}: {name}")
    assert not offenders, "macOS 10.9-incompatible libc++ calls:\n" + "\n".join(
        offenders
    )


def test_ci_builds_with_pinned_vcv_compatible_libcxx_headers():
    workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
    assert 'VCV_LIBCXX_COMMIT: "900c3b6b832d1d0e7d6e1220f6ba001802cbe0cc"' in workflow
    assert "plugin-build-mac-libcxx-compat:" in workflow
    assert "-mmacosx-version-min=10.9" in workflow
    assert '-nostdinc++ -isystem "$VCV_LIBCXX_INCLUDE"' in workflow
    assert "is unavailable: introduced in macOS 10.13" in workflow
    assert "std::any_cast<int>" in workflow
    assert "std::visit" in workflow


def test_rack_and_default_cmake_builds_do_not_require_python_analysis():
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
    cmake = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    assert "python/triggerfish_percussion" not in makefile
    assert re.search(r'option\(TRIGGERFISH_BUILD_PYTHON\s+"[^"]+"\s+OFF\)', cmake)
    assert "if(TRIGGERFISH_BUILD_PYTHON)" in cmake


def test_percussion_package_import_does_not_load_scipy():
    script = "import triggerfish_percussion, sys; assert 'scipy' not in sys.modules"
    subprocess.run([sys.executable, "-c", script], check=True, cwd=ROOT)
