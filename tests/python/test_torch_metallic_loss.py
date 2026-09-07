"""Optional analysis autograd must agree with the canonical NumPy measurements."""

import numpy as np
import pytest

torch = pytest.importorskip("torch")

from triggerfish_percussion.metallic_balance_loss import MetallicBalanceLoss
from triggerfish_percussion.torch_metallic_loss import TorchMetallicLoss


@pytest.fixture(scope="module")
def fixture():
    torch.set_num_threads(1)
    rate = 24000
    t = np.arange(rate * 4) / rate
    tone = 0.1 * np.sin(2 * np.pi * 517 * t) * np.exp(-t / 1.4)
    noise = 0.05 * np.random.default_rng(92).normal(size=len(t)) * np.exp(-t / 0.7)
    canonical = MetallicBalanceLoss(tone + noise, rate)
    return tone, noise, canonical, TorchMetallicLoss(canonical)


def test_exact_measurements(fixture):
    tone, noise, canonical, measurement = fixture
    for audio in (tone + noise, 0.5 * tone + 2 * noise):
        assert measurement.validate(audio)["relative_error"] < 1e-8


def test_gradient_matches_finite_difference(fixture):
    tone, noise, canonical, measurement = fixture
    weights = torch.tensor([0.6, 1.4], dtype=torch.float64, requires_grad=True)
    value = measurement(
        weights[0] * torch.tensor(tone) + weights[1] * torch.tensor(noise)
    )
    value.backward()
    for index in range(2):
        plus = np.array([0.6, 1.4])
        minus = plus.copy()
        plus[index] += 1e-5
        minus[index] -= 1e-5
        score = (
            lambda v: np.linalg.norm(canonical.residual(v[0] * tone + v[1] * noise))
            ** 2
        )
        numeric = (score(plus) - score(minus)) / 2e-5
        assert weights.grad[index].item() == pytest.approx(numeric, rel=1e-4, abs=1e-4)


def test_gain_optimizer_recovers_known_answer(fixture):
    from scipy.optimize import minimize

    tone, noise, _, measurement = fixture
    columns = torch.tensor(np.stack([tone, noise]))

    def objective(values):
        weights = torch.tensor(values, requires_grad=True)
        value = measurement(weights @ columns)
        value.backward()
        return value.item(), weights.grad.numpy()

    result = minimize(
        objective,
        [0.5, 1.5],
        jac=True,
        bounds=[(0.1, 2)] * 2,
        method="L-BFGS-B",
        options=dict(maxiter=20, ftol=1e-9),
    )
    np.testing.assert_allclose(result.x, [1, 1], atol=0.005)


def test_erb_weighted_measurement_parity(fixture):
    tone, noise, _, _ = fixture
    canonical = MetallicBalanceLoss(
        tone + noise, 24000, contrast_weighting="erb", fast_attack=True
    )
    assert (
        TorchMetallicLoss(canonical).validate(0.7 * tone + 1.4 * noise)[
            "relative_error"
        ]
        < 1e-8
    )


def test_silence_has_zero_loss_and_finite_gradient(fixture):
    tone, _, _, _ = fixture
    silence = np.zeros_like(tone)
    canonical = MetallicBalanceLoss(silence, 24000, fast_attack=True)
    measurement = TorchMetallicLoss(canonical)
    audio = torch.tensor(silence, requires_grad=True)
    value = measurement(audio)
    value.backward()
    # Frequency convolution can differ from SciPy by float64 roundoff.
    assert value.item() == pytest.approx(0, abs=1e-25)
    assert measurement.validate(silence)["actual"] == pytest.approx(0, abs=1e-25)
    assert torch.isfinite(audio.grad).all()
    assert torch.count_nonzero(audio.grad) == 0


def test_float32_audio_matches_canonical_and_preserves_gradient(fixture):
    tone, noise, _, _ = fixture
    reference = (tone + noise).astype(np.float32)
    candidate = (0.7 * tone + 1.4 * noise).astype(np.float32)
    canonical = MetallicBalanceLoss(reference, 24000, fast_attack=True)
    measurement = TorchMetallicLoss(canonical)
    assert measurement.validate(candidate)["relative_error"] < 1e-8
    audio = torch.tensor(candidate, requires_grad=True)
    measurement(audio).backward()
    assert audio.grad.dtype == torch.float32
    assert torch.isfinite(audio.grad).all()


@pytest.mark.parametrize("bad", [float("nan"), float("inf")])
def test_validation_fails_closed_on_nonfinite_objective(fixture, monkeypatch, bad):
    tone, noise, canonical, measurement = fixture
    monkeypatch.setattr(canonical, "residual", lambda audio: np.array([bad]))
    with pytest.raises(ValueError, match="measurement differs"):
        measurement.validate(tone + noise)


@pytest.mark.parametrize("kind", ["shape", "complex", "nan", "inf"])
def test_invalid_audio_rejected(fixture, kind):
    tone, _, _, measurement = fixture
    audio = torch.tensor(tone)
    if kind == "shape":
        audio = audio[None, :]
    elif kind == "complex":
        audio = audio.to(torch.complex128)
    else:
        audio[0] = float(kind)
    with pytest.raises(ValueError, match="finite real mono CPU"):
        measurement(audio)
