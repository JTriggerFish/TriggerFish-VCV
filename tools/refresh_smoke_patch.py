"""Prepare a private Rack smoke patch without overwriting local edits."""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path


def _device_settings(template_path: Path) -> dict[str, dict]:
    """Read only Rack MIDI/audio selections from an uncompressed local patch."""
    try:
        document = json.loads(template_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return {}
    settings = {}
    for module in document.get("modules", []):
        model = module.get("model")
        data = module.get("data", {})
        if model in {"MIDIToCVInterface", "MIDICCToCVInterface"} and isinstance(
            data.get("midi"), dict
        ):
            settings[model] = {"midi": data["midi"]}
        elif model == "AudioInterface" and isinstance(data.get("audio"), dict):
            settings[model] = {"audio": data["audio"]}
    return settings


def _portable_with_device_settings(
    portable_data: bytes, template_path: Path | None
) -> bytes:
    if template_path is None or not template_path.is_file():
        return portable_data
    settings = _device_settings(template_path)
    if not settings:
        return portable_data
    try:
        document = json.loads(portable_data.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError):
        return portable_data
    for module in document.get("modules", []):
        inherited = settings.get(module.get("model"))
        if inherited is None and module.get("model") == "MIDICCToCVInterface":
            # A refreshed electric-piano patch may be inheriting from an older
            # local fixture that predates its CC64 module. Use the keyboard
            # selected by MIDI-CV so sustain works without a second setup pass.
            inherited = settings.get("MIDIToCVInterface")
        if inherited:
            module.setdefault("data", {}).update(inherited)
    return (json.dumps(document, indent=2) + "\n").encode("utf-8")


def prepare_local_patch(
    portable_path: Path,
    local_path: Path,
    device_template: Path | None = None,
    *,
    refresh: bool = False,
) -> bool:
    """Create or explicitly refresh *local_path* from *portable_path*.

    Rack may save ``.vcv`` files as Zstandard-compressed tar archives. Existing
    local patches are treated as opaque and left byte-for-byte unchanged unless
    ``refresh`` is requested. During a refresh, readable MIDI/audio device
    selections are copied from the old local patch before its topology and
    parameters are replaced with the current portable smoke fixture.

    Returns ``True`` when an existing local patch was reused.
    """

    if local_path.exists() and not refresh:
        return True

    if refresh and device_template is None and local_path.is_file():
        device_template = local_path
    portable_data = _portable_with_device_settings(
        portable_path.read_bytes(), device_template
    )
    local_path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=local_path.parent, prefix=f".{local_path.name}.", suffix=".tmp"
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(portable_data)
        temporary_path.replace(local_path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise

    return False


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("portable_patch", type=Path)
    parser.add_argument("local_patch", type=Path)
    parser.add_argument(
        "--device-template",
        type=Path,
        help="Private local patch whose MIDI/audio selections should be inherited",
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Replace an existing local patch with the current portable fixture",
    )
    args = parser.parse_args()
    existed = args.local_patch.exists()
    reused = prepare_local_patch(
        args.portable_patch,
        args.local_patch,
        args.device_template,
        refresh=args.refresh,
    )
    action = "Using existing" if reused else "Refreshed" if existed else "Created"
    print(f"{action} private local patch {args.local_patch.name}.")


if __name__ == "__main__":
    main()
