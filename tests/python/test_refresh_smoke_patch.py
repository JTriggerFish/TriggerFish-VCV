import sys
from pathlib import Path

ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(ROOT))
from tools.refresh_smoke_patch import prepare_local_patch


def test_prepare_creates_exact_local_copy(tmp_path):
    portable_path = tmp_path / "test.vcv"
    local_path = tmp_path / "test.local.vcv"
    portable_data = b'{"modules": []}\n'
    portable_path.write_bytes(portable_data)

    assert prepare_local_patch(portable_path, local_path) is False
    assert local_path.read_bytes() == portable_data


def test_prepare_preserves_existing_compressed_local_patch(tmp_path):
    portable_path = tmp_path / "test.vcv"
    local_path = tmp_path / "test.local.vcv"
    portable_path.write_bytes(b'{"modules": []}\n')
    compressed_local_data = b"\x28\xb5\x2f\xfd\x00user Rack archive"
    local_path.write_bytes(compressed_local_data)

    assert prepare_local_patch(portable_path, local_path) is True
    assert local_path.read_bytes() == compressed_local_data
