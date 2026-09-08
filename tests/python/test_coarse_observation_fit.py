import numpy as np
import pytest

from triggerfish_percussion.coarse_observation_fit import interpolation_weights
from triggerfish_percussion.coarse_observation_fit import CoarseObservationBasis
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.spectral_bloom_basis import SpectralBloomBasis


def test_broad_shape_is_positive_partition_of_unity_and_log_frequency_linear():
    weights = interpolation_weights([10, 100, 200, 400, 800, 1600], [100, 400, 1600])
    assert np.all(weights >= 0)
    np.testing.assert_allclose(weights.sum(axis=1), 1)
    np.testing.assert_allclose(weights[2], [0.5, 0.5, 0])
    np.testing.assert_allclose(weights[4], [0, 0.5, 0.5])
    np.testing.assert_allclose(
        weights[[0, 1, 3, 5]], [[1, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
    )
    assert np.all(np.diff(weights @ [0.1, 0.5, 1]) >= 0)


@pytest.mark.parametrize("knots", [[0, 100], [100, 100], [200, 100], [100, np.nan]])
def test_bad_curve_is_rejected(knots):
    with pytest.raises(ValueError):
        interpolation_weights([100], knots)


def test_grouped_actual_render_basis_and_cached_gradient():
    class Renderer:
        def render(self, parameters, seconds, seed):
            t = np.arange(round(seconds * 8000)) / 8000
            return sum(
                10 ** (parameters[f"resolved_level_{i}"] / 20)
                * np.sin(2 * np.pi * f * t)
                * np.exp(-t * (1 + i / 8))
                for i, f in enumerate([100, 200, 400, 800, 1600, 2400])
            )

    parameters = {f"resolved_level_{i}": -72 for i in range(32)}
    for i, f in enumerate([100, 200, 400, 800, 1600, 2400]):
        parameters.update({f"resolved_level_{i}": -20, f"resolved_frequency_{i}": f})
    renderer = Renderer()
    basis = CoarseObservationBasis(
        renderer, parameters, 6, (None,), [100, 400, 1600, 3000]
    )
    amplitudes = np.array([0.12, 0.3, 0.7, 0.9])
    target = renderer.render(basis.parameters(amplitudes), 6, None)
    loss = SpectralBloomLoss(target, 8000)
    cache = SpectralBloomBasis(basis, loss)
    assert cache.validate(basis, amplitudes) < 1e-8
    probe = amplitudes * 0.8
    residual, jacobian = cache.evaluate(probe)
    step = 1e-6
    direction = np.array([0.3, -0.2, 0.4, 0.1])
    numeric = (
        cache.evaluate(probe + step * direction)[0]
        - cache.evaluate(probe - step * direction)[0]
    ) / (2 * step)
    np.testing.assert_allclose(jacobian @ direction, numeric, rtol=1e-5, atol=1e-6)
    assert np.linalg.norm(residual) > 0
    with pytest.raises(ValueError):
        basis.parameters([0, 1, 1, 1])
    irregular = dict(parameters, resolved_level_1=-3)
    with pytest.raises(ValueError, match="coarse observation curve"):
        CoarseObservationBasis(renderer, irregular, 6, (None,), [100, 400, 1600, 3000])
