"""Prepare a private Rack smoke patch without overwriting local edits."""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path


def prepare_local_patch(portable_path: Path, local_path: Path) -> bool:
    """Create *local_path* from *portable_path* when it does not exist.

    Rack may save ``.vcv`` files as Zstandard-compressed tar archives. Existing
    local patches are therefore treated as opaque, user-owned files and left
    byte-for-byte unchanged. This also preserves local topology, parameter
    edits, and MIDI/audio device selections.

    Returns ``True`` when an existing local patch was reused.
    """

    if local_path.exists():
        return True

    portable_data = portable_path.read_bytes()
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
    args = parser.parse_args()
    reused = prepare_local_patch(args.portable_patch, args.local_patch)
    action = "Using existing" if reused else "Created"
    print(f"{action} private local patch {args.local_patch.name}.")


if __name__ == "__main__":
    main()
