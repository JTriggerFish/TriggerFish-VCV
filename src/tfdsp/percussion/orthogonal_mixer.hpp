#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace tfdsp::percussion {

// In-place butterfly of Givens rotations. Every setting is orthogonal; zero
// coupling is exactly the identity, including for non-power-of-two line counts.
template <std::size_t Size> class OrthogonalMixer {
  static_assert(Size > 0, "an orthogonal mixer needs at least one channel");

public:
  using Frame = std::array<float, Size>;

  void SetCoupling(float coupling) noexcept {
    coupling = std::clamp(std::isfinite(coupling) ? coupling : 0.f, 0.f, 1.f);
    constexpr float QuarterPi = .78539816339744830962f;
    cosine_ = std::cos(QuarterPi * coupling);
    sine_ = std::sin(QuarterPi * coupling);
  }

  Frame Process(Frame frame) const noexcept {
    for (std::size_t stride = 1; stride < Size; stride *= 2) {
      for (std::size_t block = 0; block < Size; block += 2 * stride) {
        for (std::size_t offset = 0;
             offset < stride && block + offset + stride < Size; ++offset)
          Rotate(frame[block + offset], frame[block + offset + stride]);
      }
    }
    return frame;
  }

private:
  void Rotate(float &first, float &second) const noexcept {
    const float rotatedFirst = cosine_ * first + sine_ * second;
    second = -sine_ * first + cosine_ * second;
    first = rotatedFirst;
  }

  float cosine_{1.f};
  float sine_{};
};

} // namespace tfdsp::percussion
