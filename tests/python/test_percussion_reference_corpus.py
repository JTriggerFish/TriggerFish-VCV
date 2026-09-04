import importlib.util
import sys
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "percussion_reference_corpus",
    Path(__file__).parents[2] / "tools/percussion_reference_corpus.py",
)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
build_catalog = MODULE.build_catalog
reference_path = MODULE.reference_path

TOOLS = Path(__file__).parents[2] / "tools"
sys.path.insert(0, str(TOOLS))
from serve_percussion_workbench import ReferenceCatalog  # noqa: E402


def _touch(root: Path, relative: str) -> None:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"RIFF")


def test_curated_reference_catalog_is_small_and_allow_listed(tmp_path: Path) -> None:
    private = tmp_path / "references/private/cells"
    private.mkdir(parents=True)
    library = tmp_path / "library"
    references = tmp_path / "references"
    for layer in ("A", "B", "C"):
        _touch(
            library,
            f"Acoustic Drums/Snares/Yamaha Mapple CA/" f"Snare YMCA 1 {layer} 01.wav",
        )
    for layer in (1, 4, 7, 10):
        _touch(
            library,
            f"Acoustic Drums/Kicks/Yamaha Oak Custom/"
            f"Kick Yamaha Oak Custom {layer} 01.wav",
        )
    for take in (1, 2, 3, 5, 7):
        _touch(library, f"Percussion/Gongs/Gong Dresden {take:02d}.wav")
    for token in ("Cl", "H-Cl", "H-Op", "Op", "RidCrsh"):
        for layer in range(1, 5):
            _touch(
                library,
                "Acoustic Drums/Cymbals/HiHats/14 K Custom Hi-Def/"
                f"Hi-Hat 14 K Custom {token} {layer:02d}.wav",
            )
    for location in ("bell", "normal", "shoulder"):
        for dynamic in ("pp", "mf", "ff"):
            _touch(
                references,
                f"21ride.stick.{location}.{dynamic}.stereo.wav",
            )

    corpora, paths = build_catalog(private, library, references)
    counts = {corpus["id"]: len(corpus["cells"]) for corpus in corpora}
    assert counts == {
        "acoustic-snare-maple": 3,
        "acoustic-kick-oak": 4,
        "gong-dresden": 5,
        "ride-21-reference": 9,
        "hihat-14-reference": 20,
    }
    calibrations = {
        corpus["calibration"]["id"]: corpus["calibration"]
        for corpus in corpora
        if corpus.get("calibration")
    }
    assert set(calibrations) == {
        "snare-standard",
        "kick-standard",
        "gong-standard",
        "ride-standard",
        "hihat-standard",
    }
    assert calibrations["snare-standard"] == {
        "id": "snare-standard",
        "name": "Snare — medium standard hit (unverified start)",
        "recipe": "drum.snare.v1",
        "parameter_preset": "snare-default",
        "articulation": "main",
        "velocity": 82,
        "repeat": 1,
    }
    assert calibrations["kick-standard"]["recipe"] == "drum.membrane.v1"
    assert calibrations["kick-standard"]["parameter_preset"] == "acoustic-kick"
    gong = next(corpus for corpus in corpora if corpus["id"] == "gong-dresden")
    assert all(cell["implement"] == 0.5 for cell in gong["cells"])
    assert all(cell["contactSpread"] == 0.3 for cell in gong["cells"])
    hihat = next(corpus for corpus in corpora if corpus["id"] == "hihat-14-reference")
    constraints = {cell["articulation"]: cell["constraint"] for cell in hihat["cells"]}
    assert constraints == {
        "closed": 1.0,
        "half-closed": 0.72,
        "half-open": 0.38,
        "open": 0.0,
        "edge": 0.0,
    }
    for corpus in corpora:
        for cell in corpus["cells"]:
            assert reference_path(paths, cell["url"]).is_file()
    assert reference_path(paths, "/reference/acoustic-kick-oak/../secret.wav") is None
    assert all("C:" not in str(corpus) for corpus in corpora)


def test_served_reference_catalog_refreshes_atomically(tmp_path: Path) -> None:
    private = tmp_path / "references/private/cells"
    private.mkdir(parents=True)
    library = tmp_path / "library"
    references = tmp_path / "references"
    catalog = ReferenceCatalog(private, library, references)
    assert catalog.refresh() == []

    relative = "Percussion/Gongs/Gong Dresden 01.wav"
    _touch(library, relative)
    corpora = catalog.refresh()
    assert [corpus["id"] for corpus in corpora] == ["gong-dresden"]
    url = corpora[0]["cells"][0]["url"]
    assert catalog.resolve(url) == (library / relative).resolve()
