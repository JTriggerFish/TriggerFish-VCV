#pragma once

#include "tfdsp/sampleRate.hpp"

namespace tfdsp::percussion {

template <typename ResamplerType> struct FixedResamplerFactory;

template <> struct FixedResamplerFactory<tfdsp::DummyResampler> {
  static auto Create() { return tfdsp::CreateDummyResampler(); }
};

template <> struct FixedResamplerFactory<tfdsp::X2Resampler_Order7> {
  static auto Create() { return tfdsp::CreateX2Resampler_Chebychev7(); }
};

template <> struct FixedResamplerFactory<tfdsp::X4Resampler_Order7> {
  static auto Create() { return tfdsp::CreateX4Resampler_Cheby7(); }
};

} // namespace tfdsp::percussion
