"""Refresh a private Rack smoke patch while retaining local device choices."""

from __future__ import annotations

import argparse
import copy
import json
import os
import tempfile
from pathlib import Path

PRESERVED_DATA_KEYS = {
    ("Core", "MIDIToCVInterface"): ("midi",),
    ("Core", "AudioInterface"): ("audio",),
}


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _module_key(module: dict) -> tuple[object, object, object]:
    return module.get("id"), module.get("plugin"), module.get("model")


def _find_previous_module(
    module: dict, candidates: list[dict], used: set[int]
) -> dict | None:
    for index, candidate in enumerate(candidates):
        if index not in used and _module_key(candidate) == _module_key(module):
            used.add(index)
            return candidate

    module_type = module.get("plugin"), module.get("model")
    for index, candidate in enumerate(candidates):
        if (
            index not in used
            and (
                candidate.get("plugin"),
                candidate.get("model"),
            )
            == module_type
        ):
            used.add(index)
            return candidate
    return None


def refresh_local_patch(portable_path: Path, local_path: Path) -> bool:
    """Write the current portable patch to *local_path* and retain device data.

    Returns ``True`` when an existing local patch supplied device state.
    Other serialized settings deliberately come from the portable patch so its
    topology and test defaults cannot silently become stale.
    """

    portable = _load(portable_path)
    reused_device_state = local_path.exists()

    if reused_device_state:
        local = _load(local_path)
        local_modules = local.get("modules", [])
        used_local_modules: set[int] = set()
        for module in portable.get("modules", []):
            preserved_keys = PRESERVED_DATA_KEYS.get(
                (module.get("plugin"), module.get("model"))
            )
            if not preserved_keys:
                continue
            previous = _find_previous_module(module, local_modules, used_local_modules)
            if previous is None:
                continue

            previous_data = previous.get("data", {})
            current_data = module.setdefault("data", {})
            for key in preserved_keys:
                if key in previous_data:
                    current_data[key] = copy.deepcopy(previous_data[key])

    serialized = json.dumps(portable, indent=2) + "\n"
    local_path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=local_path.parent, prefix=f".{local_path.name}.", suffix=".tmp"
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as output:
            output.write(serialized)
        temporary_path.replace(local_path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise

    return reused_device_state


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("portable_patch", type=Path)
    parser.add_argument("local_patch", type=Path)
    args = parser.parse_args()
    reused = refresh_local_patch(args.portable_patch, args.local_patch)
    action = "Refreshed" if reused else "Created"
    print(f"{action} private local patch {args.local_patch.name}.")


if __name__ == "__main__":
    main()
