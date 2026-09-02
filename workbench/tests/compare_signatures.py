"""Compare native MinGW and WebAssembly renders using stable summaries."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess


def signature(command: list[str]) -> dict[str, float]:
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    return json.loads(completed.stdout.strip().splitlines()[-1])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--node", type=Path, required=True)
    parser.add_argument("--wasm-test", type=Path, required=True)
    parser.add_argument("--module", type=Path, required=True)
    parser.add_argument("--relative-tolerance", type=float, default=5.0e-4)
    arguments = parser.parse_args()

    native = signature([str(arguments.native)])
    wasm = signature(
        [str(arguments.node), str(arguments.wasm_test), str(arguments.module)]
    )
    if native["api"] != wasm["api"] or native["frames"] != wasm["frames"]:
        raise SystemExit("native and WebAssembly render contracts differ")
    for name in ("peak", "energy", "absoluteSum", "earlyEnergy"):
        scale = max(abs(native[name]), abs(wasm[name]), 1.0e-20)
        relative = abs(native[name] - wasm[name]) / scale
        if relative > arguments.relative_tolerance:
            raise SystemExit(
                f"{name} differs by {relative:.6g}; "
                f"native={native[name]:.17g}, wasm={wasm[name]:.17g}"
            )
    print(json.dumps({"native": native, "wasm": wasm}, indent=2))


if __name__ == "__main__":
    main()
