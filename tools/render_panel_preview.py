"""Render a Rack panel with its component-library controls, without running Rack.

The editable panel SVG remains the visual source. Widget positions are parsed
from the module's C++ constructor, and Rack's installed component SVGs are
embedded as data URLs in a standalone preview. A Chromium-family browser is
used for the optional PNG render when one is available.
"""

from __future__ import annotations

import argparse
import ast
import base64
import binascii
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import tempfile
import time
import xml.etree.ElementTree as ET
import zlib

ROOT = Path(__file__).resolve().parents[1]
SVG = "{http://www.w3.org/2000/svg}"

CONTROL_PATTERN = re.compile(
    r"add(?P<kind>Param|Input|Output)\s*\(\s*"
    r"create(?:LightParam|Param|Input|Output)<(?P<type>[^>]+)>\s*\(\s*"
    r"Vec\(\s*(?P<x>[0-9.]+)\s*,\s*(?P<y>[0-9.]+)\s*\)\s*,\s*module\s*,\s*"
    r"Tf303VoiceCore::(?P<id>[A-Z0-9_]+)",
    re.MULTILINE,
)

MODULE_NAMES = (
    "TfSlop",
    "TfSlop4",
    "TfVDPO",
    "TfVCA",
    "Tf303Oscillator",
    "Tf303VoiceCore",
    "Tf4072VoiceCore",
    "TfWavefoldOscillator",
    "TfUnisonOscillator",
    "TfScenePack4",
    "TfReverb",
    "TfTransport",
    "TfProgSequencer",
)


def control_pattern(module_name: str) -> re.Pattern[str]:
    return re.compile(
        r"^[ \t]*add(?P<kind>Param|Input|Output)\s*\(\s*"
        r"create(?:LightParam|Param|Input|Output)<(?P<type>[^>]+)>\s*\(\s*"
        r"Vec\(\s*(?P<x>[^,]+?)\s*,\s*(?P<y>[^)]+?)\s*\)\s*,\s*module\s*,\s*"
        + re.escape(module_name)
        + r"::(?P<id>[A-Z0-9_]+)",
        re.MULTILINE,
    )


def light_pattern(module_name: str) -> re.Pattern[str]:
    return re.compile(
        r"^[ \t]*addChild\s*\(\s*createLight<(?P<type>.+?)>\s*\(\s*"
        r"Vec\(\s*(?P<x>[^,]+?)\s*,\s*(?P<y>[^)]+?)\s*\)\s*,\s*module\s*,\s*"
        + re.escape(module_name)
        + r"::(?P<id>[A-Z0-9_]+)",
        re.MULTILINE,
    )


SCREW_PATTERN = re.compile(
    r"createWidget<ScrewSilver>\s*\(\s*Vec\(\s*(?P<x>[^,]+?)\s*,\s*"
    r"(?P<y>[^)]+?)\s*\)\s*\)",
    re.MULTILINE,
)


SWITCH_DEFAULT_PATTERN = re.compile(
    r"configSwitch\(\s*(?P<id>[A-Z0-9_]+)\s*,\s*[^,]+,\s*[^,]+,\s*"
    r"(?P<default>-?[0-9.]+)f?\s*,",
    re.MULTILINE,
)


def switch_defaults(source: str) -> dict[str, float]:
    return {
        match.group("id"): float(match.group("default"))
        for match in SWITCH_DEFAULT_PATTERN.finditer(source)
    }


COMPONENTS = {
    "TfLargeAudioKnob": (
        "Davies1900hLargeBlack_bg.svg",
        "Davies1900hLargeBlack.svg",
    ),
    "TfAudioKob": ("Davies1900hBlack_bg.svg", "Davies1900hBlack.svg"),
    "TfCvKnob": ("RoundBlackKnob_bg.svg", "RoundBlackKnob.svg"),
    "TfSnapKnob": ("RoundBlackKnob_bg.svg", "RoundBlackKnob.svg"),
    "TfRotarySwitchKnob": (
        "RoundBigBlackKnob_bg.svg",
        "RoundBigBlackKnob.svg",
    ),
    "TfTrimpot": ("Trimpot_bg.svg", "Trimpot.svg"),
    "TfSlider": ("TfSlider.svg", "TfSliderHandle.svg"),
    "TfEnvelopeSlider": ("TfSlider.svg", "TfSliderHandle.svg"),
    "CKSS": ("CKSS_0.svg", "CKSS_1.svg"),
    "CKSSThree": ("CKSSThree_0.svg", "CKSSThree_1.svg", "CKSSThree_2.svg"),
    "LEDButton": ("VCVButton_0.svg",),
    "PJ301MPort": ("PJ301M.svg",),
}

LOCAL_COMPONENT_ASSETS = {"TfSlider.svg", "TfSliderHandle.svg"}

# Representative envelope positions make the two banks legible without
# pretending that the static preview is a saved patch.
SLIDER_VALUES = {
    "FILTER_ATTACK": 0.16,
    "FILTER_DECAY": 0.46,
    "FILTER_SUSTAIN": 0.64,
    "FILTER_RELEASE": 0.40,
    "AMP_ATTACK": 0.10,
    "AMP_DECAY": 0.40,
    "AMP_SUSTAIN": 0.82,
    "AMP_RELEASE": 0.34,
    "MORPH_ALIVE": 0.38,
    "FOLD_ALIVE": 0.62,
    "SYMMETRY_ALIVE": 0.48,
}

SLIDER_LIGHT_BRIGHTNESS = {
    "FILTER_DECAY": 1.0,
    "AMP_SUSTAIN": 1.0,
}

# Rack draws the LED colour and bezel procedurally, then places the component
# SVG's reflections over it. Dimensions are Rack units (3 mm at 96 DPI).
LIGHTS = {
    "MediumLight<BlueLight>": ("MediumLight.svg", 3.0 * 96.0 / 25.4),
    "MediumLight<GreenLight>": ("MediumLight.svg", 3.0 * 96.0 / 25.4),
    "MediumLight<YellowLight>": ("MediumLight.svg", 3.0 * 96.0 / 25.4),
    "MediumLight<RedLight>": ("MediumLight.svg", 3.0 * 96.0 / 25.4),
}

PANEL_GRAPHICS = (
    (
        ROOT / "res" / "TfDiodeConnectedTransistor.svg",
        64.0,
        266.0,
        112.0,
        112.0,
        0.32,
    ),
)
MODULE_GRAPHICS = {
    "TfSlop": (),
    "TfSlop4": (),
    "TfVDPO": (),
    "TfVCA": (),
    "Tf303VoiceCore": PANEL_GRAPHICS,
    "Tf303Oscillator": ((ROOT / "res" / "logo.svg", 16.0, 232.0, 148.0, 80.8, 0.12),),
    "Tf4072VoiceCore": (),
    "TfWavefoldOscillator": (),
    "TfUnisonOscillator": (),
    "TfScenePack4": (),
    "TfReverb": (),
    "TfTransport": (),
    "TfProgSequencer": (),
}

MODULE_PANEL_FILES = {
    # Prog Sequencer is dynamically resizable; document its 30 HP default.
    "TfProgSequencer": "TfProgSequencer-30.svg",
}


def module_preview_markup(module_name: str) -> str:
    if module_name == "TfTransport":
        # TfTempoDisplay is drawn by NanoVG at runtime. Mirror its factory
        # value and geometry without baking a second label into the panel SVG.
        return (
            '<text id="transport-tempo-preview" x="60" y="110.5" '
            'text-anchor="middle" dominant-baseline="middle" '
            'font-family="Share Tech Mono, monospace" font-size="9" '
            'fill="#ffd714">120.00 BPM</text>\n'
        )
    if module_name == "TfReverb":
        # TfRoomPlanWidget is drawn by NanoVG at runtime rather than by the
        # panel SVG. Reproduce its factory appearance here so documentation
        # previews do not show an empty hole where the interactive plan lives.
        left, top, width, height = 22.0, 32.0, 196.0, 82.0

        def point(x: float, y: float) -> tuple[float, float]:
            return left + x * width, top + y * height

        grid = []
        for division in range(1, 4):
            fraction = division / 4.0
            grid.append(
                f'<path d="M{left + fraction * width:g} {top:g}v{height:g}'
                f'M{left:g} {top + fraction * height:g}h{width:g}"/>'
            )
        sources = []
        for index, position in enumerate(
            ((0.30, 0.35), (0.433333, 0.35), (0.566667, 0.35), (0.70, 0.35)),
            start=1,
        ):
            x, y = point(*position)
            sources.append(
                f'<circle cx="{x:g}" cy="{y:g}" r="5.2" fill="#ffb032" '
                'stroke="#181818" stroke-width="1"/>'
                f'<text x="{x:g}" y="{y + 2:g}" text-anchor="middle" '
                'font-family="DejaVu Sans, sans-serif" font-size="5.8" '
                f'fill="#181818">{index}</text>'
            )
        listener_x, listener_y = point(0.50, 0.682)
        return (
            '<g id="reverb-room-plan-preview" data-preview-source-count="4">'
            f'<rect x="{left:g}" y="{top:g}" width="{width:g}" '
            f'height="{height:g}" rx="2" fill="#1b2024" '
            'stroke="#101010" stroke-width="1.5"/>'
            '<g fill="none" stroke="#ffffff" stroke-opacity="0.094" '
            f'stroke-width="0.7">{"".join(grid)}</g>'
            f'<text x="{left + 5:g}" y="{top + 9:g}" '
            'font-family="DejaVu Sans, sans-serif" font-size="5.2" '
            'fill="#ffb032">SOURCES</text>'
            f'<text x="{left + width - 5:g}" y="{top + 9:g}" '
            'text-anchor="end" font-family="DejaVu Sans, sans-serif" '
            'font-size="5.2" fill="#36c8eb">LISTENER</text>'
            f'{"".join(sources)}'
            f'<circle cx="{listener_x:g}" cy="{listener_y:g}" r="4.6" '
            'fill="#36c8eb"/>'
            f'<circle cx="{listener_x:g}" cy="{listener_y:g}" r="2" '
            'fill="#1b2024"/>'
            "</g>\n"
        )
    if module_name != "TfProgSequencer":
        return ""
    lines = (
        "riff = sequence {",
        "  subdiv 8n",
        "  tonic C@4",
        "  scale minor",
        "  notes 1 x2 3{quiet} _ [>4 ^5{stacc}] 6*3 ~ 8{ten}",
        "  offset -6ms 0 +6ms |> rate 1/2",
        "  cv1 0 5 0 |> interp smooth",
        "}",
        "",
        "play riff",
    )
    text = "".join(
        f'<text x="11" y="{40 + 20 * index}" fill="#fec274">'
        f'{line.replace(">", "&gt;")}</text>'
        for index, line in enumerate(lines)
    )
    return (
        '<g id="prog-sequencer-editor-preview" '
        'font-family="Share Tech Mono, monospace" font-size="12">'
        '<rect x="5" y="23" width="385" height="334" fill="#000004"/>'
        '<rect x="5" y="359" width="385" height="16" fill="#000004"/>'
        f"{text}"
        '<text x="9" y="371" fill="#fc8c62" font-size="10">PLAY 3.00  ACTIVE</text>'
        '<circle cx="404" cy="121" r="3.1" fill="#16351f" stroke="#0b160e"/>'
        '<circle cx="404" cy="121" r="1.9" fill="#54e878"/>'
        "</g>\n"
    )


def _coordinate_variables(source: str, position: int) -> dict[str, float]:
    assignments = re.compile(
        r"(?:\b(?:auto|constexpr\s+float|float)\s+)?"
        r"(?P<name>leftMargin|spacing|offset)\s*=\s*"
        r"(?P<value>[+-]?[0-9.]+)f?\s*;"
    )
    values: dict[str, float] = {}
    for assignment in assignments.finditer(source, 0, position):
        values[assignment.group("name")] = float(assignment.group("value"))
    return values


def _evaluate_coordinate(expression: str, variables: dict[str, float]) -> float:
    expression = re.sub(r"(?<=[0-9.])f\b", "", expression.strip())
    tree = ast.parse(expression, mode="eval")

    def evaluate(node: ast.AST) -> float:
        if isinstance(node, ast.Expression):
            return evaluate(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return float(node.value)
        if isinstance(node, ast.Name) and node.id in variables:
            return variables[node.id]
        if isinstance(node, ast.BinOp) and isinstance(
            node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div)
        ):
            left = evaluate(node.left)
            right = evaluate(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            return left / right
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            value = evaluate(node.operand)
            return value if isinstance(node.op, ast.UAdd) else -value
        raise ValueError(f"Unsupported widget coordinate: {expression}")

    return evaluate(tree)


def panel_coordinate(
    expression: str,
    panel_width: float,
    panel_height: float,
    variables: dict[str, float] | None = None,
) -> float:
    expression = expression.replace("box.size.x", str(panel_width))
    expression = expression.replace("RACK_GRID_WIDTH", "15")
    expression = expression.replace("RACK_GRID_HEIGHT", str(panel_height))
    return _evaluate_coordinate(expression, variables or {})


def screw_positions(
    widget_source: str, panel_width: float, panel_height: float
) -> tuple[tuple[float, float], ...]:
    return tuple(
        (
            panel_coordinate(match.group("x"), panel_width, panel_height),
            panel_coordinate(match.group("y"), panel_width, panel_height),
        )
        for match in SCREW_PATTERN.finditer(widget_source)
    )


def centered_component_specs(
    module_name: str, widget_source: str, panel_width: float
) -> tuple[tuple[str, float, float], ...]:
    if module_name != "TfProgSequencer":
        return ()

    layout = widget_source.split("void applyPanelWidth", 1)[1].split("if (parent", 1)[0]
    variables: dict[str, float] = {}
    assignments = re.compile(
        r"const float (?P<name>right|leftColumn|rightColumn)\s*=\s*"
        r"(?P<expression>[^;]+);"
    )
    for assignment in assignments.finditer(layout):
        variables[assignment.group("name")] = panel_coordinate(
            assignment.group("expression"), panel_width, 380.0, variables
        )

    positions: list[tuple[str, float, float]] = []
    for array_name in ("inputPositions", "outputPositions"):
        array = re.search(
            rf"const Vec {array_name}\[\]\s*=\s*\{{(?P<body>.*?)\}};",
            layout,
            re.DOTALL,
        )
        if array is None:
            raise RuntimeError(f"No {array_name} layout was found in {module_name}.cpp")
        for item in re.finditer(
            r"Vec\(\s*(?P<x>[^,]+?)\s*,\s*(?P<y>[^)]+?)\s*\)",
            array.group("body"),
        ):
            positions.append(
                (
                    "PJ301MPort",
                    _evaluate_coordinate(item.group("x"), variables),
                    _evaluate_coordinate(item.group("y"), variables),
                )
            )

    expected = len(
        re.findall(r"create(?:Input|Output)Centered<PJ301MPort>", widget_source)
    )
    if len(positions) != expected:
        raise RuntimeError(
            f"{module_name} declares {expected} centred ports but its layout "
            f"contains {len(positions)} positions"
        )
    return tuple(positions)


def control_coordinates(
    control: re.Match[str], widget_source: str
) -> tuple[float, float]:
    variables = _coordinate_variables(widget_source, control.start())
    return (
        _evaluate_coordinate(control.group("x"), variables),
        _evaluate_coordinate(control.group("y"), variables),
    )


# Representative positions make the preview read naturally. They affect only
# the rendered knob indicators; parameter behavior remains defined in C++.
KNOB_ANGLES = {
    "CUTOFF": 0.0,
    "RESONANCE": -145.0,
    "DRIVE": 35.0,
    "ACCENT_SWEEP_MODE": 48.0,
    "BASS": -145.0,
    "ENV_AMOUNT": -50.0,
    "NORMAL_DECAY": -15.0,
    "ACCENT_AMOUNT": 0.0,
    "ACCENT_DECAY": -35.0,
    "VCA_DECAY": 0.0,
    "CV_AMOUNT": 0.0,
    "FM_AMOUNT": -145.0,
    "RES_AMOUNT": -145.0,
    "VCA_CV_AMOUNT": 145.0,
    "OCTAVE": 0.0,
    "TUNE": 0.0,
    "SLIDE_TIME": -25.0,
    "SHAPE": 0.0,
    "WAVE": -145.0,
    "TIME_AMOUNT": 0.0,
    "SHAPE_AMOUNT": 0.0,
    "WAVE_AMOUNT": 0.0,
    "FM_MODE": 0.0,
    "ENVELOPE_CURVE": 0.0,
}

MODULE_KNOB_ANGLES = {
    "Tf303Oscillator": {
        "FM_AMOUNT": 0.0,
        "SLIDE_TIME": 43.0,
    },
    "TfWavefoldOscillator": {
        "MORPH": 0.0,
        "FOLD": -29.0,
        "SYMMETRY": -145.0,
        "UNISON_VOICES": -145.0,
        "UNISON_SPREAD": 0.0,
    },
    "TfUnisonOscillator": {
        "PULSE_WIDTH": 0.0,
        "PULSE_WIDTH_CV_AMOUNT": 0.0,
        "PWM_RATE": -23.0,
        "SPREAD": -54.0,
        "WIDTH": 44.0,
        "SUB_LEVEL": -145.0,
        "DRIFT_SPEED": 0.0,
        "HUM": -116.0,
        "COMMON_DRIFT": -128.0,
        "INDIVIDUAL_DRIFT": -128.0,
        "TRACKING": 145.0,
        "SPREAD_CV_AMOUNT": 0.0,
        "WIDTH_CV_AMOUNT": 0.0,
    },
    "TfReverb": {
        "SPACE": 29.0,
        "ASPECT": 0.0,
        "PRE_DELAY": -131.1,
        "DECAY": 40.6,
        "DAMPING": -63.8,
        "DIFFUSION": 92.8,
        "MODULATION": -58.0,
        "SHIMMER": -145.0,
        "WIDTH": 32.8,
        "BALANCE": 0.0,
        "LOW_CUT": -145.0,
        "HIGH_CUT": 95.5,
        "MIX": -43.5,
        "LEVEL": 118.6,
    },
    "TfTransport": {
        "TEMPO": -48.3,
    },
}


def svg_dimensions(path: Path) -> tuple[float, float]:
    root = ET.parse(path).getroot()

    def pixels(value: str | None) -> float:
        if value is None:
            raise ValueError(f"SVG has no dimensions: {path}")
        return float(re.sub(r"[^0-9.+-]", "", value))

    width = root.get("width")
    height = root.get("height")
    if width is not None and height is not None and "%" not in width + height:
        return pixels(width), pixels(height)

    view_box = root.get("viewBox")
    if view_box:
        _, _, view_width, view_height = (float(value) for value in view_box.split())
        return view_width, view_height
    raise ValueError(f"SVG has no dimensions: {path}")


def embedded_image(path: Path) -> str:
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/svg+xml;base64,{encoded}"


def component_asset(component_directory: Path, asset_name: str) -> Path:
    if asset_name in LOCAL_COMPONENT_ASSETS:
        return ROOT / "res" / asset_name
    return component_directory / asset_name


def image_element(
    asset: Path,
    x: float,
    y: float,
    *,
    angle: float | None = None,
    width: float | None = None,
    height: float | None = None,
    opacity: float = 1.0,
) -> str:
    source_width, source_height = svg_dimensions(asset)
    width = source_width if width is None else width
    height = source_height if height is None else height
    transform = ""
    if angle is not None:
        transform = (
            f' transform="rotate({angle:g} {x + width / 2:g} {y + height / 2:g})"'
        )
    return (
        f'  <image x="{x:g}" y="{y:g}" width="{width:g}" height="{height:g}"'
        f' opacity="{opacity:g}"{transform} href="{embedded_image(asset)}"/>\n'
    )


def light_element(asset: Path, x: float, y: float, size: float) -> str:
    centre_x = x + size / 2.0
    centre_y = y + size / 2.0
    radius = size / 2.0
    return (
        f'  <circle cx="{centre_x:g}" cy="{centre_y:g}" r="{radius:g}" '
        'fill="#333" stroke="#000" stroke-opacity="0.21" stroke-width="1"/>\n'
        + image_element(asset, x, y, width=size, height=size)
    )


def graphic_element(spec: tuple[object, ...]) -> str:
    asset, x, y, *appearance = spec
    if not appearance:
        return image_element(asset, x, y)
    width, height, opacity = appearance
    return image_element(asset, x, y, width=width, height=height, opacity=opacity)


def find_browser() -> Path | None:
    configured = os.environ.get("PANEL_PREVIEW_BROWSER")
    candidates = [
        configured,
        shutil.which("msedge"),
        shutil.which("chrome"),
        shutil.which("chromium"),
        r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return Path(candidate)
    return None


def crop_png_height(path: Path, target_height: int) -> None:
    """Crop Chromium's unscaled RGB/RGBA screenshot without image dependencies."""
    signature = b"\x89PNG\r\n\x1a\n"
    source = path.read_bytes()
    if not source.startswith(signature):
        raise ValueError(f"Browser screenshot is not a PNG: {path}")

    chunks: list[tuple[bytes, bytes]] = []
    position = len(signature)
    while position < len(source):
        length = struct.unpack(">I", source[position : position + 4])[0]
        chunk_type = source[position + 4 : position + 8]
        chunk_data = source[position + 8 : position + 8 + length]
        chunks.append((chunk_type, chunk_data))
        position += 12 + length

    ihdr = next(data for chunk_type, data in chunks if chunk_type == b"IHDR")
    width, height, bit_depth, colour_type, compression, filtering, interlace = (
        struct.unpack(">IIBBBBB", ihdr)
    )
    if height <= target_height:
        return
    channels = {0: 1, 2: 3, 4: 2, 6: 4}.get(colour_type)
    if bit_depth != 8 or channels is None or interlace != 0:
        raise ValueError("Unsupported Chromium screenshot PNG format")

    compressed = b"".join(data for kind, data in chunks if kind == b"IDAT")
    scanlines = zlib.decompress(compressed)
    stride = 1 + width * channels
    cropped = scanlines[: stride * target_height]
    replacement_ihdr = struct.pack(
        ">IIBBBBB",
        width,
        target_height,
        bit_depth,
        colour_type,
        compression,
        filtering,
        interlace,
    )

    def png_chunk(kind: bytes, data: bytes) -> bytes:
        return (
            struct.pack(">I", len(data))
            + kind
            + data
            + struct.pack(">I", binascii.crc32(kind + data) & 0xFFFFFFFF)
        )

    output = bytearray(signature)
    inserted_pixels = False
    for kind, data in chunks:
        if kind == b"IHDR":
            output.extend(png_chunk(kind, replacement_ihdr))
        elif kind == b"IDAT":
            if not inserted_pixels:
                output.extend(png_chunk(kind, zlib.compress(cropped, level=9)))
                inserted_pixels = True
        else:
            output.extend(png_chunk(kind, data))
    path.write_bytes(output)


def render_png(
    browser: Path, preview: Path, output: Path, width: float, height: float
) -> None:
    output.unlink(missing_ok=True)
    failures = []
    # Edge occasionally refuses a second headless profile immediately after a
    # previous process exits. Retrying with a fresh profile makes --all renders
    # deterministic without sharing browser state between modules.
    for attempt in range(3):
        with tempfile.TemporaryDirectory(
            prefix="browser-", dir=output.parent, ignore_cleanup_errors=True
        ) as profile:
            command = [
                str(browser),
                "--headless",
                "--disable-gpu",
                "--hide-scrollbars",
                "--no-first-run",
                "--force-device-scale-factor=1",
                # Chromium reserves a small vertical strip in some headless
                # builds. The extra viewport avoids clipping the panel; the
                # dependency-free PNG crop below removes the blank tail.
                f"--window-size={int(2 * width)},{int(2 * height) + 64}",
                f"--user-data-dir={profile}",
                f"--screenshot={output.resolve()}",
                preview.resolve().as_uri(),
            ]
            completed = subprocess.run(
                command, check=False, capture_output=True, text=True
            )
            if completed.returncode == 0:
                deadline = time.monotonic() + 2.0
                while not output.exists() and time.monotonic() < deadline:
                    time.sleep(0.05)
                if output.exists():
                    # On Windows the Edge launcher can return before the
                    # renderer has replaced its initial screenshot file.
                    time.sleep(0.75)
        if completed.returncode == 0 and output.exists():
            crop_png_height(output, int(2 * height))
            return
        detail = completed.stderr.strip() or completed.stdout.strip()
        failures.append(f"exit {completed.returncode}: {detail or 'no output'}")
        time.sleep(0.2 * (attempt + 1))
    raise RuntimeError("browser PNG render failed; " + "; ".join(failures))


def render_preview(
    rack_runtime: Path,
    output_directory: Path,
    png: bool,
    module_name: str = "Tf303VoiceCore",
    include_components: bool = True,
) -> Path:
    if module_name not in MODULE_GRAPHICS:
        raise ValueError(f"Unsupported panel module: {module_name}")
    editable_panel = (
        ROOT / "res-src" / MODULE_PANEL_FILES.get(module_name, f"{module_name}.svg")
    )
    panel_source = (
        editable_panel
        if editable_panel.is_file()
        else ROOT / "res" / f"{module_name}.svg"
    )
    widget_source = ROOT / "src" / f"{module_name}.cpp"
    component_directory = rack_runtime / "res" / "ComponentLibrary"
    if not component_directory.is_dir():
        raise FileNotFoundError(
            f"Rack component library not found at {component_directory}. "
            "Set RACK_RUNTIME_DIR to the Rack installation directory."
        )

    panel_width, panel_height = svg_dimensions(panel_source)
    panel = panel_source.read_text(encoding="utf-8")
    panel = re.sub(
        r'width="[^"]+" height="[^"]+"',
        f'width="{2 * panel_width:g}" height="{2 * panel_height:g}"',
        panel,
        count=1,
    )
    if not re.search(r"<svg\b[^>]*\bviewBox=", panel):
        panel = panel.replace(
            "<svg ",
            f'<svg viewBox="0 0 {panel_width:g} {panel_height:g}" ',
            1,
        )
    widgets = widget_source.read_text(encoding="utf-8")
    controls = list(control_pattern(module_name).finditer(widgets))
    lights = list(light_pattern(module_name).finditer(widgets))
    static_components = centered_component_specs(module_name, widgets, panel_width)
    switch_values = switch_defaults(widgets)
    if not controls and not static_components:
        raise RuntimeError(f"No panel controls were found in {module_name}.cpp")

    overlays = ['\n<g id="module-graphics-preview">\n']
    for graphic in MODULE_GRAPHICS[module_name]:
        overlays.append(graphic_element(graphic))
    overlays.append(module_preview_markup(module_name))
    overlays.append("</g>\n")

    if include_components:
        overlays.append('\n<g id="rack-component-preview">\n')
        screw = component_directory / "ScrewSilver.svg"
        for x, y in screw_positions(widgets, panel_width, panel_height):
            overlays.append(image_element(screw, x, y))

        for control in controls:
            component_type = control.group("type")
            asset_names = COMPONENTS.get(component_type)
            if asset_names is None:
                raise KeyError(f"No preview assets configured for {component_type}")
            x, y = control_coordinates(control, widgets)
            parameter_id = control.group("id")
            if component_type in {"CKSS", "CKSSThree"}:
                state = int(round(switch_values.get(parameter_id, 0.0)))
                state = max(0, min(state, len(asset_names) - 1))
                overlays.append(
                    image_element(
                        component_asset(component_directory, asset_names[state]), x, y
                    )
                )
                continue
            if component_type in {"TfSlider", "TfEnvelopeSlider"}:
                background = component_asset(component_directory, asset_names[0])
                handle = component_asset(component_directory, asset_names[1])
                background_width, background_height = svg_dimensions(background)
                handle_width, handle_height = svg_dimensions(handle)
                value = SLIDER_VALUES.get(parameter_id, 0.5)
                handle_x = x + (background_width - handle_width) / 2.0
                handle_y = y + (1.0 - value) * (background_height - handle_height)
                overlays.append(image_element(background, x, y))
                overlays.append(image_element(handle, handle_x, handle_y))
                brightness = SLIDER_LIGHT_BRIGHTNESS.get(parameter_id, 0.0)
                if component_type == "TfEnvelopeSlider" and brightness > 0.0:
                    light_x = handle_x + (handle_width - 5.0) / 2.0
                    light_y = handle_y + (handle_height - 7.0) / 2.0
                    overlays.append(
                        f'  <rect x="{light_x - 2.0:g}" y="{light_y - 2.0:g}" '
                        f'width="9" height="11" rx="2" fill="#29b2ef" '
                        f'opacity="{0.14 * brightness:g}"/>\n'
                    )
                    overlays.append(
                        f'  <rect x="{light_x:g}" y="{light_y:g}" width="5" '
                        f'height="7" rx="0.6" fill="#29b2ef" '
                        f'opacity="{brightness:g}"/>\n'
                    )
                continue
            for index, asset_name in enumerate(asset_names):
                module_angles = MODULE_KNOB_ANGLES.get(module_name, {})
                angle = (
                    module_angles.get(parameter_id, KNOB_ANGLES.get(parameter_id))
                    if index == len(asset_names) - 1
                    else None
                )
                overlays.append(
                    image_element(
                        component_asset(component_directory, asset_name),
                        x,
                        y,
                        angle=angle,
                    )
                )
        for component_type, center_x, center_y in static_components:
            asset_names = COMPONENTS[component_type]
            for asset_name in asset_names:
                asset = component_asset(component_directory, asset_name)
                width, height = svg_dimensions(asset)
                overlays.append(
                    image_element(
                        asset, center_x - width / 2.0, center_y - height / 2.0
                    )
                )
        for light in lights:
            light_type = light.group("type")
            light_spec = LIGHTS.get(light_type)
            if light_spec is None:
                raise KeyError(f"No preview assets configured for {light_type}")
            asset_name, size = light_spec
            x, y = control_coordinates(light, widgets)
            overlays.append(light_element(component_directory / asset_name, x, y, size))
        overlays.append("</g>\n")

    output_directory.mkdir(parents=True, exist_ok=True)
    suffix = "preview" if include_components else "panel"
    preview = output_directory / f"{module_name}-{suffix}.svg"
    preview.write_text(
        panel.replace("</svg>", "".join(overlays) + "</svg>"), encoding="utf-8"
    )

    if png:
        browser = find_browser()
        if browser is None:
            raise FileNotFoundError(
                "No Chromium-family browser found. Set PANEL_PREVIEW_BROWSER or use --no-png."
            )
        render_png(
            browser,
            preview,
            output_directory / f"{module_name}-{suffix}.png",
            panel_width,
            panel_height,
        )
    return preview


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rack-runtime",
        type=Path,
        default=Path(os.environ.get("RACK_RUNTIME_DIR", ROOT.parent / "Rack2")),
    )
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=ROOT / "build" / "panel-preview",
    )
    parser.add_argument("--no-png", action="store_true")
    parser.add_argument("--panel-only", action="store_true")
    parser.add_argument("--all", action="store_true", help="Render every module")
    parser.add_argument(
        "--documentation-directory",
        type=Path,
        help="Copy rendered PNGs here using stable documentation filenames",
    )
    parser.add_argument(
        "--module",
        choices=sorted(MODULE_GRAPHICS),
        default="Tf303VoiceCore",
    )
    arguments = parser.parse_args()
    module_names = MODULE_NAMES if arguments.all else (arguments.module,)
    for module_name in module_names:
        preview = render_preview(
            arguments.rack_runtime,
            arguments.output_directory,
            not arguments.no_png,
            module_name,
            not arguments.panel_only,
        )
        print(preview)
        if arguments.documentation_directory and not arguments.no_png:
            arguments.documentation_directory.mkdir(parents=True, exist_ok=True)
            source = arguments.output_directory / f"{module_name}-preview.png"
            destination = arguments.documentation_directory / f"{module_name}.png"
            shutil.copyfile(source, destination)


if __name__ == "__main__":
    main()
