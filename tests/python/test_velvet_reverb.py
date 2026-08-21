import importlib.util
from pathlib import Path
import re

import pytest

torch = pytest.importorskip("torch")
import _triggerfish_dsp as dsp

from triggerfish_reverb.objectives import resonance_loss, resonance_per_control
from triggerfish_reverb.velvet import (
    LINE_COUNT,
    MAIN_DELAY_RATIO,
    PERMUTATIONS,
    REFERENCE_ROOM_DIMENSIONS_METRES,
    SIGNS,
    VELVET_DELAY_MS,
    VELVET_STAGE_COUNT,
    DifferentiableVelvetReverb,
    VelvetControls,
    mean_free_time_seconds,
    room_dimensions_metres,
)


def _numbers_between(text: str, first: str, last: str) -> list[float]:
    section = text.split(first, 1)[1].split(last, 1)[0]
    return [
        float(value.rstrip("f")) for value in re.findall(r"-?\d+(?:\.\d+)?f?", section)
    ]


def test_reference_has_the_current_two_stage_topology_only():
    assert LINE_COUNT == 16
    assert VELVET_STAGE_COUNT == 2
    model = DifferentiableVelvetReverb(sample_rate=4_000.0)
    assert model.base_velvet_ms.shape == (2, 16)
    assert model.permutations.shape == (3, 16)
    source = Path("python/triggerfish_reverb/velvet.py").read_text(encoding="utf-8")
    optimizer = Path("tools/optimize_velvet_reverb.py").read_text(encoding="utf-8")
    for retired in ("householder", "hybrid", "three-stage"):
        assert retired not in source.lower()
        assert retired not in optimizer.lower()


def test_python_baseline_coefficients_match_the_cpp_header():
    header = Path("src/tfdsp/late_reverb_coefficients.hpp").read_text(encoding="utf-8")
    assert "VelvetStageCount = 2" in header
    main = _numbers_between(header, "MainDelayRatio{", "};")
    velvet = _numbers_between(header, "VelvetDelayMs{{", "}};")
    permutations = _numbers_between(header, "TransformPermutation{{", "}};")
    signs = _numbers_between(header, "TransformSign{{", "}};")
    assert main == pytest.approx(MAIN_DELAY_RATIO)
    assert velvet == pytest.approx(
        [value for stage in VELVET_DELAY_MS for value in stage]
    )
    assert permutations == [value for transform in PERMUTATIONS for value in transform]
    assert signs == [value for transform in SIGNS for value in transform]


@pytest.mark.parametrize(
    "controls",
    (
        VelvetControls(space=0.0, aspect=0.5, decay=0.0, damping=0.0, diffusion=0.0),
        VelvetControls(space=0.5, aspect=0.0, decay=0.12, damping=0.45, diffusion=0.75),
        VelvetControls(space=1.0, aspect=1.0, decay=0.12, damping=1.0, diffusion=1.0),
    ),
)
def test_cpp_runtime_wall_impulse_matches_pytorch_static_response(controls):
    sample_rate = 4_000.0
    sample_count = 16_384
    impulse = dsp.late_reverb_wall_impulse(
        sample_count,
        sample_rate=sample_rate,
        space=controls.space,
        aspect=controls.aspect,
        decay=controls.decay,
        damping=controls.damping,
        diffusion=controls.diffusion,
    )
    actual = torch.fft.rfft(torch.from_numpy(impulse), dim=0)
    all_frequencies = torch.fft.rfftfreq(sample_count, 1.0 / sample_rate)
    requested = torch.tensor((20.0, 67.0, 173.0, 511.0, 997.0, 1_499.0, 1_850.0))
    bins = torch.round(requested * sample_count / sample_rate).to(torch.long)
    frequencies = all_frequencies[bins]

    model = DifferentiableVelvetReverb(sample_rate=sample_rate)
    expected, _ = model.response(frequencies, (controls,))
    torch.testing.assert_close(actual[bins], expected[0], rtol=4.0e-3, atol=8.0e-4)


def test_exporter_requires_current_architecture_and_separate_namespace():
    spec = importlib.util.spec_from_file_location(
        "export_reverb_coefficients", "tools/export_reverb_coefficients.py"
    )
    assert spec is not None and spec.loader is not None
    exporter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(exporter)
    artifact = {
        "format": exporter.ARTIFACT_FORMAT,
        "architecture": {"line_count": 16, "velvet_stage_count": 2},
        "frequency_seconds": 32.0,
        "best_step": 0,
        "best_loss": 1.0,
        "main_delay_ratio": list(MAIN_DELAY_RATIO),
        "velvet_delay_ms": [list(stage) for stage in VELVET_DELAY_MS],
        "permutations": [list(row) for row in PERMUTATIONS],
        "signs": [list(row) for row in SIGNS],
        "history": [],
    }
    exporter.validate_artifact(artifact)
    header = exporter.render_header(artifact, "0" * 64)
    assert "namespace tfdsp::late_reverb_optimized_coefficients" in header
    assert "VelvetStageCount = 2" in header
    assert "namespace tfdsp::late_reverb_coefficients" not in header
    artifact["architecture"]["velvet_stage_count"] = 3
    with pytest.raises(ValueError):
        exporter.validate_artifact(artifact)


def test_fixed_signed_hadamard_transforms_are_orthogonal():
    model = DifferentiableVelvetReverb(sample_rate=4_000.0)
    identity = torch.eye(16)
    for transform in range(3):
        matrix = model.transform(transform)
        assert torch.allclose(matrix.T @ matrix, identity, atol=2.0e-6)


def test_wall_projection_is_an_energy_preserving_embedding():
    model = DifferentiableVelvetReverb(sample_rate=4_000.0)
    identity = torch.eye(6)
    assert torch.allclose(
        model.wall_projection.T @ model.wall_projection, identity, atol=1.0e-7
    )


def test_room_geometry_and_mean_free_time_match_the_cpp_control_law():
    space = torch.tensor([0.0, 0.5, 1.0])
    aspect = torch.full_like(space, 0.5)
    dimensions = room_dimensions_metres(space, aspect)
    assert torch.allclose(dimensions[0], torch.tensor((2.8, 3.5, 2.4)))
    assert torch.allclose(dimensions[-1], torch.tensor((18.0, 25.0, 8.0)))
    assert torch.allclose(
        dimensions[1],
        torch.tensor(REFERENCE_ROOM_DIMENSIONS_METRES),
        rtol=2.0e-6,
        atol=2.0e-6,
    )
    mean_time = mean_free_time_seconds(dimensions)
    assert torch.all(torch.diff(mean_time) > 0.0)


def test_diffusion_and_room_scale_set_the_exact_integer_velvet_taps():
    model = DifferentiableVelvetReverb(sample_rate=48_000.0)
    controls = (
        VelvetControls(space=0.5, diffusion=0.0),
        VelvetControls(space=0.5, diffusion=0.75),
        VelvetControls(space=0.5, diffusion=1.0),
        VelvetControls(space=1.0, diffusion=1.0),
    )
    actual = model.velvet_delay_samples(controls)
    _, room_scale, diffusion_scale = model.room_state(controls)
    expected = torch.round(
        model.base_velvet_ms.reshape(1, 2, 16)
        * room_scale.reshape(-1, 1, 1)
        * diffusion_scale.reshape(-1, 1, 1)
        * 48.0
    ).clamp_min(1.0)
    assert torch.equal(actual, expected)
    assert actual[0].max() < actual[1].max() < actual[2].max()
    assert actual[2].max() < actual[3].max()


def test_lossless_two_stage_vfm_is_paraunitary_at_every_static_setting():
    model = DifferentiableVelvetReverb(sample_rate=4_000.0)
    frequencies = torch.tensor((31.0, 237.0, 997.0, 1_731.0))
    controls = tuple(
        VelvetControls(space=space, aspect=aspect, diffusion=diffusion)
        for space, aspect, diffusion in (
            (0.0, 0.5, 0.0),
            (0.5, 0.0, 0.75),
            (1.0, 1.0, 1.0),
        )
    )
    response = model.feedback_response(frequencies, controls, include_losses=False)
    identity = torch.eye(16, dtype=response.dtype)
    product = response.mH @ response
    assert torch.allclose(product, identity, atol=3.0e-6)


def test_zero_damping_multiband_filter_reconstructs_one_scalar_gain():
    model = DifferentiableVelvetReverb(sample_rate=48_000.0)
    frequencies = torch.tensor((20.0, 220.0, 3_500.0, 15_000.0))
    controls = (VelvetControls(decay=0.55, damping=0.0),)
    paths = torch.linspace(0.002, 0.05, 16).reshape(1, 16)
    response = model.multiband_loss_response(frequencies, paths, controls)
    decay_seconds = 0.25 * 32.0**0.55
    expected = torch.pow(torch.tensor(10.0), -3.0 * paths / decay_seconds)
    assert torch.allclose(
        response, expected.unsqueeze(1).to(response.dtype), atol=2.0e-6
    )


def test_runtime_lagrange_delay_is_exact_for_integer_taps():
    model = DifferentiableVelvetReverb(sample_rate=4_000.0)
    frequencies = torch.tensor((20.0, 500.0, 1_500.0))
    delays = torch.tensor((100.0, 127.0)).reshape(1, 2)
    response = model._fractional_delay_response(frequencies, delays)
    omega = 2.0 * torch.pi * frequencies / model.sample_rate
    expected = torch.exp(-1j * omega.reshape(-1, 1) * delays.reshape(1, -1)).unsqueeze(
        0
    )
    assert torch.allclose(response, expected, atol=2.0e-6)


def test_static_current_response_is_finite_and_differentiable():
    model = DifferentiableVelvetReverb(sample_rate=4_000.0)
    frequencies = torch.linspace(20.0, 1_900.0, 384)
    controls = (
        VelvetControls(space=0.0, aspect=0.5, decay=1.0, damping=0.0, diffusion=0.0),
        VelvetControls(space=1.0, aspect=1.0, decay=1.0, damping=1.0, diffusion=1.0),
    )
    response, internal = model.response(frequencies, controls)
    objective = response.abs().square().mean() + internal.abs().mean()
    objective.backward()
    assert response.shape == (2, frequencies.numel(), 6, 6)
    assert internal.shape == (2, frequencies.numel(), 16, 16)
    assert torch.isfinite(response).all()
    assert model.raw_main_ratios.grad is not None
    assert torch.isfinite(model.raw_main_ratios.grad).all()
    assert model.raw_main_ratios.grad.abs().sum() > 0.0
    assert float(model.main_delay_ratios().mean().detach()) == pytest.approx(1.0)


def test_resonance_loss_uses_physical_frequency_scales_not_bin_counts():
    def score(step_hz: float) -> float:
        frequencies = torch.arange(20.0, 2_000.0, step_hz)
        magnitude_db = 5.0 * torch.sin(2.0 * torch.pi * frequencies / 16.0)
        amplitude = torch.pow(10.0, magnitude_db / 20.0)
        response = amplitude.reshape(1, -1, 1, 1).expand(-1, -1, 6, 6)
        loss, _ = resonance_loss(response.to(torch.complex64), frequencies)
        return float(loss)

    coarse = score(0.25)
    fine = score(0.0625)
    assert fine == pytest.approx(coarse, rel=0.03)


def test_control_batching_has_the_exact_full_objective_gradient():
    frequencies = torch.linspace(20.0, 1_900.0, 192)
    controls = tuple(
        VelvetControls(
            space=space,
            aspect=aspect,
            decay=0.8,
            damping=damping,
            diffusion=diffusion,
        )
        for space, aspect, damping, diffusion in (
            (0.0, 0.5, 0.0, 0.0),
            (0.5, 0.0, 0.18, 0.75),
            (0.5, 1.0, 1.0, 0.75),
            (1.0, 0.5, 1.0, 1.0),
        )
    )

    def combined_response(model, selected_controls):
        wall, resolvent = model.response(frequencies, selected_controls)
        return torch.cat((wall.flatten(2), resolvent.flatten(2)), dim=2)

    full = DifferentiableVelvetReverb(sample_rate=4_000.0)
    full_loss, _ = resonance_loss(combined_response(full, controls), frequencies)
    full_loss.backward()
    expected = full.raw_main_ratios.grad.clone()

    batched = DifferentiableVelvetReverb(sample_rate=4_000.0)
    with torch.no_grad():
        reference = torch.cat(
            [
                resonance_per_control(
                    combined_response(batched, controls[first : first + 2]),
                    frequencies,
                )[0]
                for first in range(0, len(controls), 2)
            ]
        )
        weights = torch.full_like(
            reference, 1.0 / reference.numel()
        ) + 0.75 * torch.softmax(12.0 * reference, dim=0)
    for first in range(0, len(controls), 2):
        per_control, _ = resonance_per_control(
            combined_response(batched, controls[first : first + 2]), frequencies
        )
        (weights[first : first + 2] * per_control).sum().backward()
    assert torch.allclose(
        batched.raw_main_ratios.grad, expected, rtol=3.0e-4, atol=3.0e-5
    )
