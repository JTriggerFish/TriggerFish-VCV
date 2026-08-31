#include "tfdsp/late_reverb.hpp"
#include "tfdsp/multiband_decay_filter.hpp"
#include "tfdsp/windowed_pitch_shifter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double ReverbTestSampleRate = 48'000.0;

void Check(const bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(EXIT_FAILURE);
  }
}

void TestDecayFilterRejectsInvalidT60() {
  tfdsp::MultibandDecayFilter filter;
  filter.Prepare(48000.0);
  const float infinity = std::numeric_limits<float>::infinity();
  Check(std::abs(filter.Process(.75f, .01f, infinity, infinity, infinity) -
                 .75f) < 1.e-6f,
        "positive-infinite T60 remains the explicit lossless setting");
  filter.Reset();
  const float nan = std::numeric_limits<float>::quiet_NaN();
  Check(filter.Process(.75f, .01f, nan, nan, nan) == 0.f,
        "NaN T60 cannot silently enable lossless feedback");
  filter.Reset();
  Check(filter.Process(.75f, .01f, -infinity, -infinity, -infinity) == 0.f,
        "negative-infinite T60 cannot silently enable lossless feedback");
  filter.Reset();
  Check(filter.Process(std::numeric_limits<float>::denorm_min(), .01f,
                       1.f, 1.f, 1.f) == 0.f,
        "reverb decay filter flushes subnormal state to exact silence");
}

void TestCompleteLatePathsFlushSubnormals() {
  const float denormal = std::numeric_limits<float>::denorm_min();
  tfdsp::WindowedPitchShifter shifter;
  shifter.Prepare(48000.0);
  bool shifterSilent = true;
  for (std::size_t sample = 0; sample < 8192; ++sample)
    shifterSilent = shifterSilent &&
        shifter.Process(sample == 0 ? denormal : 0.f) == 0.f;
  Check(shifterSilent,
        "windowed pitch shifter flushes subnormal filter and delay state");

  tfdsp::LateReverb reverb;
  tfdsp::LateReverbControls controls;
  controls.shimmer = 1.f;
  bool reverbSilent = true;
  for (std::size_t sample = 0; sample < 8192; ++sample) {
    const auto output = reverb.Process(sample == 0 ? denormal : 0.f, controls);
    reverbSilent = reverbSilent && output[0] == 0.f && output[1] == 0.f;
  }
  Check(reverbSilent,
        "complete late-reverb feedback and shimmer paths flush subnormals");
}

void TestLateReverbAtEverySupportedRate() {
  for (const double sampleRate : {44'100.0, 48'000.0, 88'200.0,
                                  96'000.0, 192'000.0}) {
    tfdsp::LateReverb reverb;
    reverb.SetSampleRate(sampleRate);
    tfdsp::LateReverbControls controls;
    bool finite = true;
    const auto count = static_cast<std::size_t>(sampleRate * .08);
    for (std::size_t sample = 0; sample < count; ++sample) {
      const auto output = reverb.Process(sample == 0 ? 1.f : 0.f, controls);
      finite = finite && std::isfinite(output[0]) && std::isfinite(output[1]);
    }
    Check(finite, "late reverb remains finite at every supported sample rate");
  }
}

std::array<float, 3> RoomDimensions(const float space) {
  const float x = std::clamp(space, 0.f, 1.f);
  const float amount = x * x * (3.f - 2.f * x);
  const auto interpolate = [amount](const float minimum, const float maximum) {
    return std::exp(std::log(minimum) +
                    amount * (std::log(maximum) - std::log(minimum)));
  };
  return {interpolate(2.8f, 18.f), interpolate(3.5f, 25.f),
          interpolate(2.4f, 8.f)};
}

std::vector<tfdsp::LateReverb::StereoFrame>
Render(const tfdsp::LateReverbControls &controls, const std::size_t samples,
       const double sampleRate = ReverbTestSampleRate,
       const tfdsp::LateReverbFlavour flavour =
           tfdsp::DefaultLateReverbFlavour) {
  tfdsp::LateReverb reverb;
  reverb.SetFlavour(flavour);
  reverb.SetSampleRate(sampleRate);
  std::vector<tfdsp::LateReverb::StereoFrame> result(samples);
  for (std::size_t sample = 0; sample < samples; ++sample)
    result[sample] = reverb.Process(sample == 0 ? 1.f : 0.f, controls);
  return result;
}

double Energy(const std::vector<tfdsp::LateReverb::StereoFrame> &response,
              const std::size_t first) {
  double energy = 0.0;
  for (std::size_t sample = first; sample < response.size(); ++sample)
    for (const float value : response[sample])
      energy += static_cast<double>(value) * value;
  return energy;
}

double
SegmentEnergy(const std::vector<tfdsp::LateReverb::StereoFrame> &response,
              const std::size_t first, const std::size_t last) {
  double energy = 0.0;
  for (std::size_t sample = first; sample < std::min(last, response.size());
       ++sample)
    for (const float value : response[sample])
      energy += static_cast<double>(value) * value;
  return energy;
}

double
ChannelEnergy(const std::vector<tfdsp::LateReverb::StereoFrame> &response,
              const std::size_t channel) {
  double energy = 0.0;
  for (const auto &frame : response)
    energy += static_cast<double>(frame[channel]) * frame[channel];
  return energy;
}

double BandEnergy(const std::vector<tfdsp::LateReverb::StereoFrame> &response,
                  const double sampleRate, const double lowerFrequency,
                  const double upperFrequency, const std::size_t first) {
  std::array<double, 2> belowLow{};
  std::array<double, 2> belowHigh{};
  const double lowAlpha = 1.0 - std::exp(-2.0 * 3.14159265358979323846 *
                                         lowerFrequency / sampleRate);
  const double highAlpha = 1.0 - std::exp(-2.0 * 3.14159265358979323846 *
                                          upperFrequency / sampleRate);
  double energy = 0.0;
  for (std::size_t sample = 0; sample < response.size(); ++sample)
    for (std::size_t channel = 0; channel < 2; ++channel) {
      belowLow[channel] +=
          lowAlpha * (response[sample][channel] - belowLow[channel]);
      belowHigh[channel] +=
          highAlpha * (response[sample][channel] - belowHigh[channel]);
      if (sample >= first) {
        const double band = belowHigh[channel] - belowLow[channel];
        energy += band * band;
      }
    }
  return energy;
}

std::size_t FirstAudibleSample(
    const std::vector<tfdsp::LateReverb::StereoFrame> &response) {
  for (std::size_t sample = 0; sample < response.size(); ++sample)
    if (std::abs(response[sample][0]) + std::abs(response[sample][1]) > 1.e-7f)
      return sample;
  return response.size();
}

double MidbandT60(const std::vector<tfdsp::LateReverb::StereoFrame> &response,
                  const double sampleRate) {
  std::vector<double> energy(response.size());
  std::array<double, 2> lowState{};
  std::array<double, 2> bandState{};
  const double highpassAlpha =
      1.0 - std::exp(-2.0 * 3.14159265358979323846 * 400.0 / sampleRate);
  const double lowpassAlpha =
      1.0 - std::exp(-2.0 * 3.14159265358979323846 * 4'000.0 / sampleRate);
  for (std::size_t sample = 0; sample < response.size(); ++sample)
    for (std::size_t channel = 0; channel < 2; ++channel) {
      lowState[channel] +=
          highpassAlpha * (response[sample][channel] - lowState[channel]);
      const double highpassed = response[sample][channel] - lowState[channel];
      bandState[channel] += lowpassAlpha * (highpassed - bandState[channel]);
      energy[sample] += bandState[channel] * bandState[channel];
    }
  for (std::size_t sample = energy.size() - 1; sample > 0; --sample)
    energy[sample - 1] += energy[sample];
  const double reference = std::max(energy.front(), 1.e-30);
  const auto crossing = [&](const double decibels) {
    const double ratio = std::pow(10.0, decibels / 10.0);
    for (std::size_t sample = 0; sample < energy.size(); ++sample)
      if (energy[sample] / reference <= ratio)
        return static_cast<double>(sample) / sampleRate;
    return static_cast<double>(energy.size()) / sampleRate;
  };
  return 2.0 * (crossing(-35.0) - crossing(-5.0));
}

double LateGaussianDeviation(
    const std::vector<tfdsp::LateReverb::StereoFrame> &response,
    const std::size_t first, const std::size_t last, const std::size_t window) {
  double deviation = 0.0;
  std::size_t windows = 0;
  for (std::size_t start = first;
       start + window <= std::min(last, response.size()); start += window) {
    for (std::size_t channel = 0; channel < 2; ++channel) {
      double mean = 0.0;
      for (std::size_t sample = start; sample < start + window; ++sample)
        mean += response[sample][channel];
      mean /= static_cast<double>(window);
      double second = 0.0;
      double fourth = 0.0;
      for (std::size_t sample = start; sample < start + window; ++sample) {
        const double centred = response[sample][channel] - mean;
        const double square = centred * centred;
        second += square;
        fourth += square * square;
      }
      second /= static_cast<double>(window);
      fourth /= static_cast<double>(window);
      if (second > 1.e-30) {
        deviation += std::abs(fourth / (3.0 * second * second) - 1.0);
        ++windows;
      }
    }
  }
  return deviation / std::max<std::size_t>(windows, 1);
}

double TonePower(const std::vector<float> &signal, const std::size_t first,
                 const double frequency, const double sampleRate) {
  double real = 0.0;
  double imaginary = 0.0;
  for (std::size_t sample = first; sample < signal.size(); ++sample) {
    const double phase = 2.0 * 3.14159265358979323846 * frequency *
                         static_cast<double>(sample) / sampleRate;
    real += signal[sample] * std::cos(phase);
    imaginary -= signal[sample] * std::sin(phase);
  }
  return real * real + imaginary * imaginary;
}

double ToneBandPower(const std::vector<float> &signal, const std::size_t first,
                     const double lowFrequency, const double highFrequency,
                     const double sampleRate, const double step = 0.5) {
  double power = 0.0;
  for (double frequency = lowFrequency; frequency <= highFrequency;
       frequency += step)
    power += TonePower(signal, first, frequency, sampleRate);
  return power;
}

double SignalEnergy(const std::vector<float> &signal, const std::size_t first) {
  double energy = 0.0;
  for (std::size_t sample = first; sample < signal.size(); ++sample)
    energy += static_cast<double>(signal[sample]) * signal[sample];
  return energy;
}

double LowFrequencyPeriodicity(
    const std::vector<tfdsp::LateReverb::StereoFrame> &response,
    const double sampleRate, const double lowerFrequency = 20.0,
    const double upperFrequency = 250.0,
    const double analysisStartSeconds = 0.15) {
  std::vector<double> energy(response.size());
  std::array<double, 2> lowState{};
  std::array<double, 2> bandState{};
  const double highpassAlpha = 1.0 - std::exp(-2.0 * 3.14159265358979323846 *
                                              lowerFrequency / sampleRate);
  const double lowpassAlpha = 1.0 - std::exp(-2.0 * 3.14159265358979323846 *
                                             upperFrequency / sampleRate);
  for (std::size_t sample = 0; sample < response.size(); ++sample) {
    for (std::size_t channel = 0; channel < 2; ++channel) {
      lowState[channel] +=
          highpassAlpha * (response[sample][channel] - lowState[channel]);
      const double highpassed = response[sample][channel] - lowState[channel];
      bandState[channel] += lowpassAlpha * (highpassed - bandState[channel]);
      energy[sample] += bandState[channel] * bandState[channel];
    }
  }

  const std::size_t first =
      static_cast<std::size_t>(analysisStartSeconds * sampleRate);
  const std::size_t window = static_cast<std::size_t>(0.005 * sampleRate);
  const std::size_t hop = std::max<std::size_t>(1, window / 2);
  std::vector<double> envelope;
  for (std::size_t start = first; start + window <= energy.size();
       start += hop) {
    double sum = 0.0;
    for (std::size_t sample = start; sample < start + window; ++sample)
      sum += energy[sample];
    envelope.push_back(std::log(sum / static_cast<double>(window) + 1.e-20));
  }
  const std::size_t trendRadius =
      static_cast<std::size_t>(0.05 * sampleRate / static_cast<double>(hop));
  std::vector<double> residual(envelope.size());
  for (std::size_t index = 0; index < envelope.size(); ++index) {
    const std::size_t begin = index > trendRadius ? index - trendRadius : 0;
    const std::size_t end = std::min(envelope.size(), index + trendRadius + 1);
    double trend = 0.0;
    for (std::size_t sample = begin; sample < end; ++sample)
      trend += envelope[sample];
    residual[index] =
        envelope[index] - trend / static_cast<double>(end - begin);
  }
  double mean = 0.0;
  for (const double value : residual)
    mean += value;
  mean /= residual.size();
  for (double &value : residual)
    value -= mean;

  const std::size_t minimumLag =
      static_cast<std::size_t>(0.015 * sampleRate / static_cast<double>(hop));
  const std::size_t maximumLag =
      static_cast<std::size_t>(0.200 * sampleRate / static_cast<double>(hop));
  double maximumCorrelation = 0.0;
  for (std::size_t lag = minimumLag; lag <= maximumLag; ++lag) {
    double product = 0.0;
    double earlyEnergy = 0.0;
    double lateEnergy = 0.0;
    for (std::size_t index = 0; index + lag < residual.size(); ++index) {
      product += residual[index] * residual[index + lag];
      earlyEnergy += residual[index] * residual[index];
      lateEnergy += residual[index + lag] * residual[index + lag];
    }
    maximumCorrelation = std::max(
        maximumCorrelation,
        product / std::sqrt(std::max(earlyEnergy * lateEnergy, 1.e-20)));
  }
  return maximumCorrelation;
}

double SizeAutomationFarSidebandRatio(const float diffusion,
                                      const double frequency) {
  constexpr double sampleRate = ReverbTestSampleRate;
  constexpr std::size_t sampleCount = 192'000;
  constexpr std::size_t automationStart = 96'000;
  constexpr std::size_t analysisEnd = 144'000;
  const double smoothingCoefficient =
      1.0 - std::exp(-1.0 / (0.020 * sampleRate));
  tfdsp::LateReverb reverb;
  reverb.SetSampleRate(sampleRate);
  tfdsp::LateReverbControls controls;
  float sizeControl = 0.5f;
  controls.roomDimensionsMetres = RoomDimensions(sizeControl);
  controls.decay = 1.f;
  controls.damping = 0.18f;
  controls.diffusion = diffusion;
  controls.modulation = 0.f;
  controls.shimmer = 0.f;
  std::vector<float> output;
  output.reserve(analysisEnd - automationStart);
  for (std::size_t sample = 0; sample < sampleCount; ++sample) {
    if (sample >= automationStart)
      sizeControl +=
          static_cast<float>(smoothingCoefficient) * (1.f - sizeControl);
    controls.roomDimensionsMetres = RoomDimensions(sizeControl);
    const float input =
        static_cast<float>(0.1 * std::sin(2.0 * 3.14159265358979323846 *
                                          frequency * sample / sampleRate));
    const auto frame = reverb.Process(input, controls);
    if (sample >= automationStart && sample < analysisEnd)
      output.push_back(frame[0] + frame[1]);
  }

  // A fixed LTI reverb can change a tone's amplitude and phase but cannot
  // translate it. A moving delay read does, leaving strong energy far from the
  // excitation frequency. A Hann window keeps the measurement insensitive to
  // the finite analysis interval.
  for (std::size_t sample = 0; sample < output.size(); ++sample) {
    const double phase = 2.0 * 3.14159265358979323846 * sample /
                         static_cast<double>(output.size() - 1);
    output[sample] *= static_cast<float>(0.5 - 0.5 * std::cos(phase));
  }
  double total = 0.0;
  double far = 0.0;
  for (double binFrequency = 20.0; binFrequency <= 1'900.0;
       binFrequency += 1.0) {
    const double power = TonePower(output, 0, binFrequency, sampleRate);
    total += power;
    if (std::abs(binFrequency - frequency) > 35.0)
      far += power;
  }
  return far / std::max(total, 1.e-20);
}

void TestVelvetFeedbackMatrixIsParaunitaryAndDense() {
  constexpr double sampleRate = 48'000.0;
  for (const auto flavour : {tfdsp::LateReverbFlavour::Base,
                             tfdsp::LateReverbFlavour::Optimized}) {
    for (const float diffusion : {0.f, 0.5f, 1.f}) {
      tfdsp::VelvetFeedbackMatrix matrix;
      matrix.SetFlavour(flavour);
      matrix.Prepare(sampleRate);
      tfdsp::VelvetFeedbackMatrix::Frame input{};
      input[0] = 1.f;
      double outputEnergy = 0.0;
      std::size_t nonzero = 0;
      for (std::size_t sample = 0; sample < 4'096; ++sample) {
        const auto output = matrix.Process(input, diffusion);
        input = {};
        for (const float value : output) {
          outputEnergy += static_cast<double>(value) * value;
          nonzero += std::abs(value) > 1.e-7f;
        }
      }
      Check(std::abs(outputEnergy - 1.0) < 2.e-5,
            "both velvet feedback flavours must conserve energy; error=" +
                std::to_string(std::abs(outputEnergy - 1.0)));
      if (diffusion == 1.f)
        Check(nonzero >= 32,
              "both multi-stage velvet flavours must create many sparse "
              "paths");
    }
  }
}

void TestVelvetFractionalModulationPreservesTheStaticPath() {
  constexpr float infinity = std::numeric_limits<float>::infinity();
  for (const double sampleRate : {44'100.0, 48'000.0, 96'000.0}) {
    tfdsp::VelvetFeedbackMatrix reference;
    tfdsp::VelvetFeedbackMatrix explicitZero;
    tfdsp::VelvetFeedbackMatrix modulated;
    reference.Prepare(sampleRate);
    explicitZero.Prepare(sampleRate);
    modulated.Prepare(sampleRate);
    tfdsp::VelvetFeedbackMatrix::Frame input{};
    input[0] = 1.f;
    double zeroDifference = 0.0;
    double modulatedDifference = 0.0;
    float modulatedPeak = 0.f;
    const auto sampleCount =
        static_cast<std::size_t>(std::ceil(0.1 * sampleRate));
    for (std::size_t sample = 0; sample < sampleCount; ++sample) {
      const auto expected = reference.Process(input, 0.75f);
      const auto zero = explicitZero.Process(input, 0.75f, infinity, infinity,
                                             infinity, 1.f, 0.f);
      const auto moving = modulated.Process(input, 0.75f, infinity, infinity,
                                            infinity, 1.f, 1.f);
      input = {};
      for (std::size_t line = 0; line < expected.size(); ++line) {
        zeroDifference += std::abs(expected[line] - zero[line]);
        modulatedDifference += std::abs(expected[line] - moving[line]);
        Check(std::isfinite(moving[line]),
              "fractionally modulated velvet delays must remain finite");
        modulatedPeak = std::max(modulatedPeak, std::abs(moving[line]));
      }
    }
    Check(zeroDifference == 0.0,
          "zero velvet modulation must preserve the exact integer-delay path");
    Check(modulatedDifference > 1.e-4,
          "maximum velvet modulation must move the diffusion-stage taps");
    Check(modulatedPeak < 2.f,
          "fractionally modulated velvet delays must remain bounded");
  }
}

void TestLateReverbFlavourSelectionAndTransition() {
  constexpr double sampleRate = 48'000.0;
  Check(tfdsp::LateReverb{}.Flavour() ==
            tfdsp::LateReverbFlavour::Optimized,
        "Optimized must be the default late-tail flavour");

  tfdsp::LateReverb selectedBeforePrepare;
  selectedBeforePrepare.SetFlavour(tfdsp::LateReverbFlavour::Base);
  selectedBeforePrepare.SetSampleRate(sampleRate);
  tfdsp::LateReverb restoredAfterPrepare;
  restoredAfterPrepare.SetSampleRate(sampleRate);
  restoredAfterPrepare.SetFlavourImmediate(tfdsp::LateReverbFlavour::Base);
  tfdsp::LateReverbControls restoreControls;
  restoreControls.modulation = 0.f;
  restoreControls.shimmer = 0.f;
  for (std::size_t sample = 0; sample < 12'000; ++sample) {
    const float input = sample == 0 ? 1.f : 0.f;
    const auto expected =
        selectedBeforePrepare.Process(input, restoreControls);
    const auto restored = restoredAfterPrepare.Process(input, restoreControls);
    for (std::size_t channel = 0; channel < expected.size(); ++channel)
      Check(std::abs(expected[channel] - restored[channel]) < 1.e-7f,
            "restoring a saved Base flavour before audio must not render an "
            "Optimized-to-Base startup transition");
  }

  tfdsp::LateReverbControls controls;
  controls.roomDimensionsMetres = RoomDimensions(0.75f);
  controls.decay = 0.75f;
  controls.diffusion = 0.8f;
  controls.modulation = 0.f;
  tfdsp::LateReverb switched;
  switched.SetFlavour(tfdsp::LateReverbFlavour::Base);
  switched.SetSampleRate(sampleRate);
  double beforeEnergy = 0.0;
  double afterEnergy = 0.0;
  for (std::size_t sample = 0; sample < 192'000; ++sample) {
    if (sample == 48'000)
      switched.SetFlavour(tfdsp::LateReverbFlavour::Optimized);
    const auto frame = switched.Process(sample == 0 ? 1.f : 0.f, controls);
    for (const float value : frame) {
      Check(std::isfinite(value),
            "switching late-tail flavour must remain finite");
      if (sample >= 42'000 && sample < 48'000)
        beforeEnergy += static_cast<double>(value) * value;
      if (sample >= 48'000 && sample < 54'000)
        afterEnergy += static_cast<double>(value) * value;
    }
  }
  Check(afterEnergy > 0.01 * beforeEnergy,
        "switching flavour must preserve the already-stored tail; ratio=" +
            std::to_string(afterEnergy / std::max(beforeEnergy, 1.e-20)));
  Check(switched.Flavour() == tfdsp::LateReverbFlavour::Optimized,
        "the selected Optimized flavour must remain active");
}

void TestDiffusionControlsTemporalSpread() {
  constexpr double sampleRate = 48'000.0;
  const auto support = [=](const float diffusion) {
    tfdsp::VelvetFeedbackMatrix matrix;
    matrix.Prepare(sampleRate);
    tfdsp::VelvetFeedbackMatrix::Frame input{};
    input[0] = 1.f;
    std::size_t last = 0;
    for (std::size_t sample = 0; sample < 4'096; ++sample) {
      const auto output = matrix.Process(input, diffusion);
      input = {};
      for (const float value : output)
        if (std::abs(value) > 1.e-7f)
          last = sample;
    }
    return last;
  };
  const auto minimum = support(0.f);
  const auto maximum = support(1.f);
  Check(maximum > 4 * minimum,
        "Diffusion must expand the VFM's temporal scattering span; minimum=" +
            std::to_string(minimum) + ", maximum=" + std::to_string(maximum));
}

void TestRoomScaleControlsTheCompleteVelvetSpan() {
  constexpr double sampleRate = 48'000.0;
  const auto support = [=](const float roomScale) {
    tfdsp::VelvetFeedbackMatrix matrix;
    matrix.Prepare(sampleRate);
    tfdsp::VelvetFeedbackMatrix::Frame input{};
    input[0] = 1.f;
    std::size_t last = 0;
    double energy = 0.0;
    for (std::size_t sample = 0; sample < 8'192; ++sample) {
      const auto output = matrix.Process(input, 0.75f, roomScale);
      input = {};
      for (const float value : output) {
        energy += static_cast<double>(value) * value;
        if (std::abs(value) > 1.e-7f)
          last = sample;
      }
    }
    Check(std::abs(energy - 1.0) < 2.e-5,
          "room-scaled velvet operator must remain paraunitary");
    return last;
  };
  const auto small = support(0.45f);
  const auto nominal = support(1.f);
  const auto large = support(2.15f);
  Check(small < nominal && nominal < large && large > 3 * small,
        "Size must scale both velvet delay stages; support=" +
            std::to_string(small) + "/" + std::to_string(nominal) + "/" +
            std::to_string(large));
}

void TestRoomGeometryControlsLateArrivalTime() {
  constexpr double sampleRate = 48'000.0;
  tfdsp::LateReverbControls small;
  small.roomDimensionsMetres = RoomDimensions(0.f);
  small.modulation = 0.f;
  tfdsp::LateReverbControls large = small;
  large.roomDimensionsMetres = RoomDimensions(1.f);
  const auto smallResponse = Render(small, 12'000, sampleRate);
  const auto largeResponse = Render(large, 12'000, sampleRate);
  const auto smallArrival = FirstAudibleSample(smallResponse);
  const auto largeArrival = FirstAudibleSample(largeResponse);
  Check(largeArrival > smallArrival * 3 / 2,
        "larger room geometry must produce a later late-field arrival; small=" +
            std::to_string(smallArrival) +
            ", large=" + std::to_string(largeArrival));
}

void TestDampingProducesAudibleFrequencyDependentDecay() {
  constexpr double sampleRate = 48'000.0;
  constexpr std::size_t sampleCount = 192'000;
  tfdsp::LateReverbControls bright;
  bright.decay = 0.75f;
  bright.damping = 0.f;
  bright.modulation = 0.f;
  tfdsp::LateReverbControls damped = bright;
  damped.damping = 1.f;
  const auto brightResponse = Render(bright, sampleCount, sampleRate);
  const auto dampedResponse = Render(damped, sampleCount, sampleRate);
  constexpr std::size_t tailStart = 48'000;
  const double brightLow =
      BandEnergy(brightResponse, sampleRate, 250.0, 1'000.0, tailStart);
  const double dampedLow =
      BandEnergy(dampedResponse, sampleRate, 250.0, 1'000.0, tailStart);
  const double brightHigh =
      BandEnergy(brightResponse, sampleRate, 5'000.0, 12'000.0, tailStart);
  const double dampedHigh =
      BandEnergy(dampedResponse, sampleRate, 5'000.0, 12'000.0, tailStart);
  const double lowRatio = dampedLow / std::max(brightLow, 1.e-20);
  const double highRatio = dampedHigh / std::max(brightHigh, 1.e-20);
  Check(highRatio < 0.20 * lowRatio,
        "Damping must shorten high-frequency decay much more than low-mid "
        "decay; high ratio=" +
            std::to_string(highRatio) +
            ", low ratio=" + std::to_string(lowRatio));
}

void TestMaximumRoomLongDecayAndModulationRemainStableWithoutShimmer() {
  constexpr double sampleRate = 48'000.0;
  constexpr std::size_t sampleCount = 576'000;
  tfdsp::LateReverbControls controls;
  controls.roomDimensionsMetres = RoomDimensions(1.f);
  controls.decay = 1.f;
  controls.damping = 0.18f;
  controls.modulation = 1.f;
  controls.shimmer = 0.f;
  for (const auto flavour :
       {tfdsp::LateReverbFlavour::Base, tfdsp::LateReverbFlavour::Optimized}) {
    for (const float diffusion : {0.f, 0.75f, 1.f}) {
      controls.diffusion = diffusion;
      const auto response = Render(controls, sampleCount, sampleRate, flavour);
      float peak = 0.f;
      for (const auto &frame : response)
        for (const float value : frame) {
          Check(std::isfinite(value),
                "maximum-room, maximum-modulation tail must remain finite "
                "without shimmer at every Diffusion setting");
          peak = std::max(peak, std::abs(value));
        }
      const double middleTail = SegmentEnergy(response, 192'000, 288'000);
      const double finalTail = SegmentEnergy(response, 480'000, 576'000);
      const double gaussianDeviation =
          LateGaussianDeviation(response, 96'000, 480'000, 2'400);
      double maximumPeriodicity = 0.0;
      double problemBand = 0.0;
      for (const auto band : {std::array<double, 2>{20.0, 80.0},
                              std::array<double, 2>{80.0, 160.0},
                              std::array<double, 2>{160.0, 320.0},
                              std::array<double, 2>{320.0, 640.0},
                              std::array<double, 2>{640.0, 1'280.0},
                              std::array<double, 2>{1'280.0, 2'560.0}}) {
        const double periodicity = LowFrequencyPeriodicity(
            response, sampleRate, band[0], band[1], 2.0);
        if (periodicity > maximumPeriodicity) {
          maximumPeriodicity = periodicity;
          problemBand = band[0];
        }
      }
      Check(peak < 2.f,
            "maximum-room, maximum-modulation tail must retain safe signal "
            "headroom");
      Check(finalTail < 0.5 * middleTail,
            "maximum-room tail must decay at every Diffusion setting; ratio=" +
                std::to_string(finalTail / std::max(middleTail, 1.e-20)) +
                ", diffusion=" + std::to_string(diffusion));
      Check(maximumPeriodicity < 0.15,
            "maximum-room tail must not form repeating octave-band energy "
            "packets; periodicity=" +
                std::to_string(maximumPeriodicity) +
                ", band-start=" + std::to_string(problemBand) +
                ", diffusion=" + std::to_string(diffusion));
      Check(gaussianDeviation < 0.20,
            "maximum-room tail must retain a dense noise-like distribution; "
            "Gaussian deviation=" +
                std::to_string(gaussianDeviation) +
                ", diffusion=" + std::to_string(diffusion));
    }
  }
}

void TestImpulseIsFiniteDenseAndStereo() {
  constexpr double sampleRate = 48'000.0;
  const auto response = Render({}, 96'000, sampleRate);
  std::size_t nonzero = 0;
  bool channelsDiffer = false;
  float peak = 0.f;
  for (const auto &frame : response) {
    for (const float value : frame) {
      Check(std::isfinite(value), "late response must remain finite");
      peak = std::max(peak, std::abs(value));
      nonzero += std::abs(value) > 1.e-8f;
    }
    channelsDiffer |= std::abs(frame[0] - frame[1]) > 1.e-7f;
  }
  Check(nonzero > 60'000, "late response should become temporally dense");
  Check(channelsDiffer, "default response should be stereo");
  Check(peak > 1.e-5f && peak < 2.f,
        "late response should have bounded useful level");
  const double imbalanceDb = 10.0 * std::log10(ChannelEnergy(response, 0) /
                                               ChannelEnergy(response, 1));
  Check(std::abs(imbalanceDb) < 0.5,
        "default native late response must not be biased toward either ear; "
        "imbalance=" +
            std::to_string(imbalanceDb) + " dB");
}

void TestDefaultBassTailAvoidsPeriodicEchoBuildUp() {
  constexpr double sampleRate = 48'000.0;
  tfdsp::LateReverbControls controls;
  controls.decay = 0.55f;
  controls.damping = 0.18f;
  controls.diffusion = 0.75f;
  const auto response = Render(controls, 96'000, sampleRate);
  const double periodicity = LowFrequencyPeriodicity(response, sampleRate);
  Check(periodicity < 0.15, "default low-frequency tail should not form "
                            "periodic echo packets; score=" +
                                std::to_string(periodicity));
}

void TestDecayControl() {
  tfdsp::LateReverbControls shortControls;
  shortControls.decay = 0.f;
  tfdsp::LateReverbControls longControls = shortControls;
  longControls.decay = 1.f;
  const auto shortResponse = Render(shortControls, 144'000);
  const auto longResponse = Render(longControls, 144'000);
  Check(Energy(longResponse, 48'000) > 20.0 * Energy(shortResponse, 48'000),
        "long decay must retain substantially more late energy");
}

void TestMidbandT60TracksDecayControl() {
  constexpr double sampleRate = 48'000.0;
  for (const float decay : {0.4f, 0.6f}) {
    tfdsp::LateReverbControls controls;
    controls.roomDimensionsMetres = RoomDimensions(0.5f);
    controls.decay = decay;
    controls.damping = 0.f;
    controls.diffusion = 0.75f;
    controls.modulation = 0.f;
    controls.shimmer = 0.f;
    const double target = 0.25 * std::exp(decay * std::log(32.0));
    const auto response =
        Render(controls, static_cast<std::size_t>(4.0 * target * sampleRate),
               sampleRate);
    const double measured = MidbandT60(response, sampleRate);
    Check(std::abs(measured / target - 1.0) < 0.12,
          "fixed VFM midband T60 must follow the Decay control; target=" +
              std::to_string(target) +
              " s, measured=" + std::to_string(measured) + " s");
  }
}

void TestMidbandT60RemainsIndependentOfRoomSize() {
  constexpr double sampleRate = ReverbTestSampleRate;
  constexpr float decay = 0.55f;
  const double target = 0.25 * std::exp(decay * std::log(32.0));
  for (const float space : {0.f, 0.5f, 1.f}) {
    tfdsp::LateReverbControls controls;
    controls.roomDimensionsMetres = RoomDimensions(space);
    controls.decay = decay;
    controls.damping = 0.f;
    controls.diffusion = 0.75f;
    controls.modulation = 0.f;
    controls.shimmer = 0.f;
    const auto response =
        Render(controls, static_cast<std::size_t>(4.0 * target * sampleRate),
               sampleRate);
    const double measured = MidbandT60(response, sampleRate);
    Check(std::abs(measured / target - 1.0) < 0.15,
          "Decay must retain its RT60 while Size changes; size=" +
              std::to_string(space) + ", target=" + std::to_string(target) +
              " s, measured=" + std::to_string(measured) + " s");
  }
}

void TestWindowedPitchShifterRaisesAnOctave() {
  constexpr double sampleRate = 48'000.0;
  constexpr double inputFrequency = 440.0;
  tfdsp::WindowedPitchShifter shifter;
  shifter.Prepare(sampleRate, 0.120f);
  Check(std::abs(shifter.WindowSamples() / sampleRate - 0.120) < 1.e-5,
        "pitch-shift grain duration must be sample-rate invariant");

  std::vector<float> output(144'000);
  float peak = 0.f;
  for (std::size_t sample = 0; sample < output.size(); ++sample) {
    const float input = static_cast<float>(
        std::sin(2.0 * 3.14159265358979323846 * inputFrequency *
                 static_cast<double>(sample) / sampleRate));
    output[sample] = shifter.Process(input);
    Check(std::isfinite(output[sample]), "pitch shifter must remain finite");
    peak = std::max(peak, std::abs(output[sample]));
  }
  const std::size_t settled = 48'000;
  const double octavePower =
      TonePower(output, settled, 2.0 * inputFrequency, sampleRate);
  const double originalPower =
      TonePower(output, settled, inputFrequency, sampleRate);
  Check(octavePower > 20.0 * originalPower,
        "granular resampling must move a sustained tone up one octave");
  Check(peak < 1.5f, "pitch-shift crossfade must keep bounded amplitude");

  const double centrePower =
      ToneBandPower(output, settled, 877.0, 883.0, sampleRate);
  const double grainRate = sampleRate / shifter.WindowSamples();
  const double periodicSidebandPower =
      ToneBandPower(output, settled, 880.0 - grainRate - 2.0,
                    880.0 - grainRate + 2.0, sampleRate) +
      ToneBandPower(output, settled, 880.0 + grainRate - 2.0,
                    880.0 + grainRate + 2.0, sampleRate);
  Check(periodicSidebandPower < 0.6 * centrePower,
        "randomized grains must not concentrate more energy in the grain-rate "
        "sidebands than around the intended octave; ratio=" +
            std::to_string(periodicSidebandPower /
                           std::max(centrePower, 1.e-20)));
}

void TestPitchShifterAdvancePreservesWarmState() {
  constexpr double sampleRate = 48'000.0;
  constexpr std::size_t warmSamples = 48'000;
  constexpr std::size_t comparisonSamples = 12'000;
  tfdsp::WindowedPitchShifter rendered;
  tfdsp::WindowedPitchShifter advanced;
  rendered.Prepare(sampleRate, 0.120f, 0.23f, 0x1234abcdU);
  advanced.Prepare(sampleRate, 0.120f, 0.23f, 0x1234abcdU);

  for (std::size_t sample = 0; sample < warmSamples + comparisonSamples;
       ++sample) {
    const float input = static_cast<float>(
        0.7 * std::sin(2.0 * 3.14159265358979323846 * 317.0 * sample /
                       sampleRate) +
        0.2 * std::sin(2.0 * 3.14159265358979323846 * 1'831.0 * sample /
                       sampleRate));
    const float expected = rendered.Process(input);
    if (sample < warmSamples) {
      advanced.Advance(input);
    } else {
      const float actual = advanced.Process(input);
      Check(std::abs(actual - expected) < 1.e-7f,
            "the cheap zero-gain shimmer path must preserve complete pitch-"
            "shifter state");
    }
  }
}

void TestPitchShifterHasNoGrainBoundarySpikes() {
  constexpr double sampleRate = 48'000.0;
  constexpr double inputFrequency = 440.0;
  constexpr std::size_t settled = 48'000;
  constexpr std::size_t sampleCount = 192'000;
  tfdsp::WindowedPitchShifter shifter;
  shifter.Prepare(sampleRate, 0.120f);

  float previous = 0.f;
  double differenceEnergy = 0.0;
  float peakDifference = 0.f;
  std::size_t differenceCount = 0;
  for (std::size_t sample = 0; sample < sampleCount; ++sample) {
    const float input = static_cast<float>(
        std::sin(2.0 * 3.14159265358979323846 * inputFrequency *
                 static_cast<double>(sample) / sampleRate));
    const float output = shifter.Process(input);
    if (sample > settled) {
      const float difference = output - previous;
      differenceEnergy += static_cast<double>(difference) * difference;
      peakDifference = std::max(peakDifference, std::abs(difference));
      ++differenceCount;
    }
    previous = output;
  }

  const double differenceRms =
      std::sqrt(differenceEnergy / std::max<std::size_t>(differenceCount, 1));
  const double crest = peakDifference / std::max(differenceRms, 1.e-20);
  Check(crest < 6.0,
        "uniform grain overlap must not expose click-like derivative spikes; "
        "difference crest factor=" +
            std::to_string(crest));
}

void TestShimmerProjectionIsOrthonormal() {
  using namespace tfdsp::late_reverb_coefficients;
  for (std::size_t first = 0; first < ShimmerBusCount; ++first)
    for (std::size_t second = 0; second < ShimmerBusCount; ++second) {
      double innerProduct = 0.0;
      for (std::size_t line = 0; line < LineCount; ++line)
        innerProduct += ShimmerProjection[first][line] *
                        ShimmerProjection[second][line];
      const double expected = first == second ? 1.0 : 0.0;
      Check(std::abs(innerProduct - expected) < 1.e-6,
            "shimmer projection rows must be orthonormal so subspace "
            "routing has deterministic unit gain");
    }
}

void TestFixedStereoFdnVectorsAreOrthonormal() {
  using namespace tfdsp::late_reverb_coefficients;
  double inputNorm = 0.0;
  std::array<double, 2> outputNorm{};
  double outputInnerProduct = 0.0;
  for (std::size_t line = 0; line < LineCount; ++line) {
    inputNorm += FixedInputVector[line] * FixedInputVector[line];
    for (std::size_t channel = 0; channel < 2; ++channel)
      outputNorm[channel] +=
          FixedStereoOutput[channel][line] *
          FixedStereoOutput[channel][line];
    outputInnerProduct +=
        FixedStereoOutput[0][line] * FixedStereoOutput[1][line];
  }
  Check(std::abs(inputNorm - 1.0) < 1.e-6,
        "fixed FDN input vector must have unit energy");
  Check(std::abs(outputNorm[0] - 1.0) < 1.e-6 &&
            std::abs(outputNorm[1] - 1.0) < 1.e-6 &&
            std::abs(outputInnerProduct) < 1.e-6,
        "fixed stereo FDN decoder rows must be orthonormal");
}

void TestShimmerReturnsOnlyThroughTheMainFeedbackLoop() {
  constexpr double sampleRate = ReverbTestSampleRate;
  constexpr std::size_t sampleCount = 24'000;
  tfdsp::LateReverb withoutShimmer;
  tfdsp::LateReverb withShimmer;
  withoutShimmer.SetSampleRate(sampleRate);
  withShimmer.SetSampleRate(sampleRate);
  tfdsp::LateReverbControls controls;
  controls.modulation = 0.f;
  controls.shimmer = 0.f;
  auto shimmerControls = controls;
  shimmerControls.shimmer = 1.f;

  std::size_t firstDifference = sampleCount;
  for (std::size_t sample = 0; sample < sampleCount; ++sample) {
    tfdsp::LateReverb::WallFrame input{};
    if (sample == 0)
      input[0] = 1.f;
    const auto normal = withoutShimmer.ProcessWallFrame(input, controls);
    const auto shimmer =
        withShimmer.ProcessWallFrame(input, shimmerControls);
    double difference = 0.0;
    for (std::size_t wall = 0; wall < normal.size(); ++wall)
      difference += std::abs(shimmer[wall] - normal[wall]);
    if (difference > 1.e-8 && firstDifference == sampleCount)
      firstDifference = sample;
  }

  const auto dimensions = controls.roomDimensionsMetres;
  const double volume = dimensions[0] * dimensions[1] * dimensions[2];
  const double surface = 2.0 * (dimensions[0] * dimensions[1] +
                                dimensions[0] * dimensions[2] +
                                dimensions[1] * dimensions[2]);
  const double meanFreeSamples =
      4.0 * volume / (343.0 * surface) * sampleRate;
  const auto minimumRatio = *std::min_element(
      tfdsp::LateMainRatios(tfdsp::DefaultLateReverbFlavour).begin(),
      tfdsp::LateMainRatios(tfdsp::DefaultLateReverbFlavour).end());
  const auto minimumRecursiveReturn = static_cast<std::size_t>(
      std::floor(2.0 * meanFreeSamples * minimumRatio));
  Check(firstDifference < sampleCount,
        "maximum shimmer must alter the recursive wall response");
  Check(firstDifference + 2 >= minimumRecursiveReturn,
        "shimmer must not bypass the main room delays or appear as a direct "
        "decoder-side layer; first difference=" +
            std::to_string(firstDifference) +
            ", minimum recursive return=" +
            std::to_string(minimumRecursiveReturn));
}

void TestPitchShifterRejectsAliasingInput() {
  constexpr double sampleRate = 48'000.0;
  constexpr std::size_t sampleCount = 96'000;
  constexpr std::size_t settled = 24'000;
  const auto render = [=](const double frequency) {
    tfdsp::WindowedPitchShifter shifter;
    shifter.Prepare(sampleRate, 0.120f);
    std::vector<float> output(sampleCount);
    for (std::size_t sample = 0; sample < output.size(); ++sample) {
      const float input = static_cast<float>(
          std::sin(2.0 * 3.14159265358979323846 * frequency *
                   static_cast<double>(sample) / sampleRate));
      output[sample] = shifter.Process(input);
    }
    return SignalEnergy(output, settled);
  };

  const double legalBandEnergy = render(4'000.0);
  const double rejectedBandEnergy = render(14'000.0);
  Check(rejectedBandEnergy < 1.e-4 * legalBandEnergy,
        "octave shifter must strongly reject input that would alias; ratio=" +
            std::to_string(rejectedBandEnergy /
                           std::max(legalBandEnergy, 1.e-20)));
}

void TestShimmerAddsAnOctaveLayerWithoutMutingTheTail() {
  tfdsp::LateReverbControls dryFeedback;
  dryFeedback.decay = 0.9f;
  dryFeedback.shimmer = 0.f;
  tfdsp::LateReverbControls shimmer = dryFeedback;
  shimmer.shimmer = 1.f;
  const auto normalResponse = Render(dryFeedback, 288'000);
  const auto shimmerResponse = Render(shimmer, 288'000);
  double difference = 0.0;
  float peak = 0.f;
  for (std::size_t sample = 4'800; sample < shimmerResponse.size(); ++sample)
    for (std::size_t channel = 0; channel < 2; ++channel) {
      difference += std::abs(shimmerResponse[sample][channel] -
                             normalResponse[sample][channel]);
      peak = std::max(peak, std::abs(shimmerResponse[sample][channel]));
      Check(std::isfinite(shimmerResponse[sample][channel]),
            "maximum shimmer output must remain finite");
    }
  const double meanDifference =
      difference /
      (2.0 * static_cast<double>(shimmerResponse.size() - 4'800));
  // Preserve the original 4 kHz test's 0.1 / (23,600 frames * 2 channels)
  // per-sample requirement now that this test renders the production rate.
  Check(meanDifference > 0.1 / 47'200.0,
        "shimmer must decisively alter the late response; mean difference=" +
            std::to_string(meanDifference));
  Check(peak < 2.f, "maximum shimmer output must remain bounded");
  const double normalTailEnergy = Energy(normalResponse, 4'800);
  const double shimmerTailEnergy = Energy(shimmerResponse, 4'800);
  Check(shimmerTailEnergy >= 0.5 * normalTailEnergy,
        "shimmer must not silence the unshifted late tail; energy ratio=" +
            std::to_string(shimmerTailEnergy /
                           std::max(normalTailEnergy, 1.e-20)));
  const double imbalanceDb =
      10.0 * std::log10(ChannelEnergy(shimmerResponse, 0) /
                        ChannelEnergy(shimmerResponse, 1));
  Check(std::abs(imbalanceDb) < 1.0,
        "shimmer output must retain a balanced stereo field");
}

void TestShimmerCreatesAnOctaveInsideTheLateField() {
  constexpr double sampleRate = ReverbTestSampleRate;
  constexpr double inputFrequency = 440.0;
  constexpr std::size_t sampleCount = 192'000;
  constexpr std::size_t settled = 48'000;
  const auto renderTone = [=](const float shimmer) {
    tfdsp::LateReverb reverb;
    reverb.SetSampleRate(sampleRate);
    tfdsp::LateReverbControls controls;
    controls.decay = 0.8f;
    controls.damping = 0.1f;
    controls.modulation = 0.f;
    controls.shimmer = shimmer;
    std::vector<float> output(sampleCount);
    for (std::size_t sample = 0; sample < output.size(); ++sample) {
      const float input = static_cast<float>(
          std::sin(2.0 * 3.14159265358979323846 * inputFrequency *
                   static_cast<double>(sample) / sampleRate));
      const auto frame = reverb.Process(input, controls);
      output[sample] = frame[0] + frame[1];
    }
    return output;
  };
  std::array<double, 5> octavePower{};
  std::array<double, 5> fundamentalPower{};
  const std::array<float, 5> shimmerValues{0.f, 0.25f, 0.5f, 0.75f, 1.f};
  for (std::size_t index = 0; index < shimmerValues.size(); ++index) {
    const auto output = renderTone(shimmerValues[index]);
    octavePower[index] =
        TonePower(output, settled, 2.0 * inputFrequency, sampleRate);
    fundamentalPower[index] =
        TonePower(output, settled, inputFrequency, sampleRate);
  }
  Check(octavePower.back() > 100.0 * std::max(octavePower.front(), 1.e-20),
        "shimmer must introduce octave-up energy into the late field");
  const double maximumOctaveRatio =
      octavePower.back() / std::max(fundamentalPower.back(), 1.e-20);
  // In a complete diffuse tail the unshifted fundamental remains dominant by
  // design. At production rate, -27 dB is an audible octave layer without
  // demanding the unsafe feedback gain implied by the former 8 kHz proxy.
  Check(maximumOctaveRatio > 0.002,
        "maximum shimmer must produce a clearly audible octave layer; "
        "octave/fundamental power ratio=" +
            std::to_string(maximumOctaveRatio));
  for (std::size_t index = 1; index < octavePower.size(); ++index) {
    Check(octavePower[index] > octavePower[index - 1],
          "shimmer must monotonically increase octave energy");
    Check(fundamentalPower[index] >= 0.8 * fundamentalPower.front(),
          "shimmer must preserve the unshifted late-field fundamental");
  }
}

void TestMaximumShimmerDoesNotGrowAfterExcitation() {
  constexpr double sampleRate = ReverbTestSampleRate;
  constexpr std::size_t sampleCount = 960'000;
  constexpr std::size_t excitationSamples = 48'000;
  constexpr std::array<std::array<float, 2>, 5> Cases{{
      {{0.f, 0.f}}, {{0.f, 1.f}}, {{0.5f, 0.75f}}, {{1.f, 0.f}}, {{1.f, 1.f}}
  }};
  for (const auto &testCase : Cases) {
    tfdsp::LateReverb reverb;
    reverb.SetSampleRate(sampleRate);
    tfdsp::LateReverbControls controls;
    controls.roomDimensionsMetres = RoomDimensions(testCase[0]);
    controls.decay = 1.f;
    controls.damping = 0.f;
    controls.diffusion = testCase[1];
    controls.shimmer = 1.f;
    std::vector<tfdsp::LateReverb::StereoFrame> response(sampleCount);
    float peak = 0.f;
    for (std::size_t sample = 0; sample < response.size(); ++sample) {
      const float input =
          sample < excitationSamples
              ? static_cast<float>(
                    0.25 * std::sin(2.0 * 3.14159265358979323846 * 220.0 *
                                    static_cast<double>(sample) / sampleRate))
              : 0.f;
      response[sample] = reverb.Process(input, controls);
      for (const float value : response[sample]) {
        Check(std::isfinite(value),
              "maximum shimmer must remain finite during a long render");
        peak = std::max(peak, std::abs(value));
      }
    }

    const double middleTail = SegmentEnergy(response, 384'000, 576'000);
    const double finalTail = SegmentEnergy(response, 768'000, 960'000);
    const std::string setting =
        "; size=" + std::to_string(testCase[0]) +
        ", diffusion=" + std::to_string(testCase[1]);
    const double decayRatio = finalTail / std::max(middleTail, 1.e-20);
    Check(peak < 4.f,
          "maximum shimmer must retain safe long-render headroom; peak=" +
              std::to_string(peak) + ", final/middle=" +
              std::to_string(decayRatio) + setting);
    Check(finalTail < middleTail,
          "maximum-shimmer tail must decay rather than enter a growing "
          "oscillation; ratio=" +
              std::to_string(decayRatio) + setting);
  }
}

void TestMaximumShimmerIsStableAtProductionSampleRate() {
  constexpr double sampleRate = 48'000.0;
  constexpr std::size_t sampleCount = 384'000;
  constexpr std::size_t excitationSamples = 24'000;
  tfdsp::LateReverb reverb;
  reverb.SetSampleRate(sampleRate);
  tfdsp::LateReverbControls controls;
  controls.roomDimensionsMetres = RoomDimensions(1.f);
  controls.decay = 1.f;
  controls.damping = 0.f;
  controls.diffusion = 1.f;
  controls.modulation = 1.f;
  controls.shimmer = 1.f;

  double middleEnergy = 0.0;
  double finalEnergy = 0.0;
  float peak = 0.f;
  for (std::size_t sample = 0; sample < sampleCount; ++sample) {
    float input = 0.f;
    if (sample < excitationSamples) {
      const double time = static_cast<double>(sample) / sampleRate;
      input = static_cast<float>(
          0.08 * (std::sin(2.0 * 3.14159265358979323846 * 173.0 * time) +
                  std::sin(2.0 * 3.14159265358979323846 * 587.0 * time) +
                  std::sin(2.0 * 3.14159265358979323846 * 1'747.0 * time) +
                  std::sin(2.0 * 3.14159265358979323846 * 5'231.0 * time)));
    }
    const auto output = reverb.Process(input, controls);
    for (const float value : output) {
      Check(std::isfinite(value),
            "production-rate maximum shimmer must remain finite");
      peak = std::max(peak, std::abs(value));
      const double square = static_cast<double>(value) * value;
      if (sample >= 192'000 && sample < 288'000)
        middleEnergy += square;
      if (sample >= 288'000)
        finalEnergy += square;
    }
  }

  Check(peak < 2.f,
        "production-rate maximum shimmer must retain safe headroom");
  Check(finalEnergy < 0.5 * middleEnergy,
        "production-rate maximum shimmer must decay after broadband "
        "excitation; final/middle energy=" +
            std::to_string(finalEnergy / std::max(middleEnergy, 1.e-20)));
}

void TestShimmerLayerIsDiffuseAndNonPeriodic() {
  constexpr double sampleRate = ReverbTestSampleRate;
  constexpr std::size_t sampleCount = 288'000;
  tfdsp::LateReverbControls controls;
  controls.decay = 0.8f;
  controls.damping = 0.18f;
  controls.diffusion = 0.75f;
  controls.modulation = 0.f;
  controls.shimmer = 0.f;
  const auto dry = Render(controls, sampleCount, sampleRate);
  controls.shimmer = 1.f;
  const auto wet = Render(controls, sampleCount, sampleRate);
  auto layer = wet;
  for (std::size_t sample = 0; sample < layer.size(); ++sample)
    for (std::size_t channel = 0; channel < 2; ++channel)
      layer[sample][channel] -= dry[sample][channel];

  double worstPeriodicity = 0.0;
  for (const auto band : {std::array<double, 2>{250.0, 500.0},
                          std::array<double, 2>{500.0, 1'000.0},
                          std::array<double, 2>{1'000.0, 2'000.0},
                          std::array<double, 2>{2'000.0, 4'000.0}})
    worstPeriodicity =
        std::max(worstPeriodicity,
                 LowFrequencyPeriodicity(layer, sampleRate, band[0], band[1],
                                         0.35));
  const double gaussian =
      LateGaussianDeviation(layer, 48'000, 192'000, 2'400);
  Check(worstPeriodicity < 0.25,
        "in-loop shimmer must not form repeating spectral-energy packets; "
        "periodicity=" +
            std::to_string(worstPeriodicity));
  Check(gaussian < 0.40,
        "in-loop shimmer must retain a dense noise-like distribution; "
        "Gaussian deviation=" +
            std::to_string(gaussian));
}

void TestSizeAutomationDoesNotDopplerTheLongTail() {
  double maximumRatio = 0.0;
  float problemDiffusion = 0.f;
  double problemFrequency = 0.0;
  for (const float diffusion : {0.f, 0.75f, 1.f})
    for (const double frequency : {220.0, 440.0, 880.0}) {
      const double ratio = SizeAutomationFarSidebandRatio(diffusion, frequency);
      if (ratio > maximumRatio) {
        maximumRatio = ratio;
        problemDiffusion = diffusion;
        problemFrequency = frequency;
      }
    }
  Check(maximumRatio < 0.02,
        "turning Size up at maximum Decay must not Doppler-shift persistent "
        "tail energy, including with non-zero Diffusion; far-sideband ratio=" +
            std::to_string(maximumRatio) +
            ", diffusion=" + std::to_string(problemDiffusion) +
            ", frequency=" + std::to_string(problemFrequency));
}

void TestFixedStereoEncodingAndModulationRangeAffectResponse() {
  tfdsp::LateReverbControls baseline;
  baseline.modulation = 0.f;
  tfdsp::LateReverbControls lowModulation = baseline;
  lowModulation.modulation = 0.1f;
  tfdsp::LateReverbControls maximumModulation = baseline;
  maximumModulation.modulation = 1.f;
  tfdsp::LateReverbControls lowDiffusion = baseline;
  lowDiffusion.diffusion = 0.f;
  tfdsp::LateReverbControls highDiffusion = baseline;
  highDiffusion.diffusion = 1.f;
  const auto baselineResponse = Render(baseline, 36'000);
  const auto repeatedResponse = Render(baseline, 36'000);
  const auto lowModulationResponse = Render(lowModulation, 36'000);
  const auto maximumModulationResponse = Render(maximumModulation, 36'000);
  const auto lowDiffusionResponse = Render(lowDiffusion, 36'000);
  const auto highDiffusionResponse = Render(highDiffusion, 36'000);
  double repeatedDifference = 0.0;
  double lowModulationDifference = 0.0;
  double maximumModulationDifference = 0.0;
  double diffusionDifference = 0.0;
  for (std::size_t sample = 0; sample < baselineResponse.size(); ++sample)
    for (std::size_t channel = 0; channel < 2; ++channel) {
      repeatedDifference += std::abs(baselineResponse[sample][channel] -
                                     repeatedResponse[sample][channel]);
      lowModulationDifference +=
          std::abs(baselineResponse[sample][channel] -
                   lowModulationResponse[sample][channel]);
      maximumModulationDifference +=
          std::abs(baselineResponse[sample][channel] -
                   maximumModulationResponse[sample][channel]);
      diffusionDifference += std::abs(lowDiffusionResponse[sample][channel] -
                                      highDiffusionResponse[sample][channel]);
    }
  const double valueCount =
      2.0 * static_cast<double>(baselineResponse.size());
  repeatedDifference /= valueCount;
  lowModulationDifference /= valueCount;
  maximumModulationDifference /= valueCount;
  diffusionDifference /= valueCount;
  Check(repeatedDifference == 0.0,
        "fixed FDN encoding is deterministic and geometry independent");
  Check(lowModulationDifference > 1.e-8,
        "modulation must have no inactive region above zero");
  Check(maximumModulationDifference > lowModulationDifference,
        "maximum modulation must move the late field more than a low setting");
  Check(diffusionDifference > 0.01 / 6'000.0,
        "Diffusion must decisively change input and feedback scattering");
}

void TestAutomatedLateControlsPreserveStoredTail() {
  constexpr double sampleRate = ReverbTestSampleRate;
  constexpr std::size_t sampleCount = 72'000;
  constexpr std::size_t excitationSamples = 12'000;
  constexpr std::size_t changeSample = 36'000;
  constexpr std::size_t measurementEnd = 48'000;

  tfdsp::LateReverbControls baseline;
  baseline.decay = 0.8f;
  baseline.damping = 0.2f;
  baseline.diffusion = 0.75f;
  baseline.modulation = 0.f;
  baseline.shimmer = 0.f;
  std::array<std::pair<const char *, tfdsp::LateReverbControls>, 6> cases{};
  for (auto &entry : cases)
    entry.second = baseline;
  cases[0].first = "room dimensions";
  cases[0].second.roomDimensionsMetres = RoomDimensions(1.f);
  cases[1].first = "Decay";
  cases[1].second.decay = 0.6f;
  cases[2].first = "Damping";
  cases[2].second.damping = 0.9f;
  cases[3].first = "Diffusion";
  cases[3].second.diffusion = 0.f;
  cases[4].first = "Modulation";
  cases[4].second.modulation = 1.f;
  cases[5].first = "Shimmer";
  cases[5].second.shimmer = 1.f;

  for (const auto &[name, target] : cases) {
    tfdsp::LateReverb reference;
    tfdsp::LateReverb changed;
    reference.SetSampleRate(sampleRate);
    changed.SetSampleRate(sampleRate);
    double referenceEnergy = 0.0;
    double changedEnergy = 0.0;
    std::uint32_t noise = 0x9e3779b9u;
    for (std::size_t sample = 0; sample < sampleCount; ++sample) {
      noise ^= noise << 13;
      noise ^= noise >> 17;
      noise ^= noise << 5;
      const float input =
          sample < excitationSamples
              ? 0.25f * (2.f * static_cast<float>(noise & 0x00ffffffu) /
                             static_cast<float>(0x00ffffffu) -
                         1.f)
              : 0.f;
      const auto referenceFrame = reference.Process(input, baseline);
      const auto changedFrame =
          changed.Process(input, sample < changeSample ? baseline : target);
      if (sample >= changeSample && sample < measurementEnd)
        for (std::size_t channel = 0; channel < 2; ++channel) {
          referenceEnergy += static_cast<double>(referenceFrame[channel]) *
                             referenceFrame[channel];
          changedEnergy += static_cast<double>(changedFrame[channel]) *
                           changedFrame[channel];
        }
    }
    const double ratio = changedEnergy / std::max(referenceEnergy, 1.e-20);
    Check(ratio > 0.05,
          std::string(name) +
              " automation must preserve the already-stored late tail; "
              "energy ratio=" + std::to_string(ratio));
    Check(ratio < 20.0,
          std::string(name) +
              " automation must remain bounded while the stored tail moves; "
              "energy ratio=" + std::to_string(ratio));
  }
}

std::vector<tfdsp::LateReverb::StereoFrame> RenderSmoke303Excitation() {
  constexpr double sampleRate = 48'000.0;
  constexpr std::size_t sampleCount = 576'000;
  constexpr std::size_t excitationSamples = 384'000;
  constexpr std::size_t stepSamples = 6'000;
  constexpr std::array<int, 16> semitones{0,  0, 12, 7,  0, 3,  10, 0,
                                          12, 7, 3,  -5, 0, 10, 7,  3};
  tfdsp::LateReverb reverb;
  reverb.SetSampleRate(sampleRate);
  tfdsp::LateReverbControls controls;
  // The checked-in smoke patch mirrors the authoritative factory scene. This
  // diagnostic therefore follows those defaults instead of maintaining a
  // second set of literals that can silently drift.
  controls.roomDimensionsMetres =
      tfdsp::reverb_defaults::RoomDimensionsMetres;
  controls.decay = tfdsp::reverb_defaults::Decay;
  controls.damping = tfdsp::reverb_defaults::Damping;
  controls.diffusion = tfdsp::reverb_defaults::Diffusion;
  controls.modulation = tfdsp::reverb_defaults::Modulation;
  controls.shimmer = tfdsp::reverb_defaults::Shimmer;
  std::vector<tfdsp::LateReverb::StereoFrame> response(sampleCount);
  double phase = 0.0;
  double filter = 0.0;
  for (std::size_t sample = 0; sample < sampleCount; ++sample) {
    float input = 0.f;
    if (sample < excitationSamples) {
      const std::size_t step = sample / stepSamples;
      const std::size_t within = sample % stepSamples;
      const double frequency =
          65.4063913 * std::exp2(semitones[step % semitones.size()] / 12.0);
      phase += frequency / sampleRate;
      phase -= std::floor(phase);
      const double envelope =
          std::exp(-static_cast<double>(within) / (0.045 * sampleRate));
      const double cutoff = 250.0 + 3'500.0 * envelope;
      const double alpha =
          1.0 - std::exp(-2.0 * 3.14159265358979323846 * cutoff / sampleRate);
      filter += alpha * ((2.0 * phase - 1.0) - filter);
      input = static_cast<float>(1.5 * envelope * filter);
    }
    response[sample] = reverb.Process(input, controls);
  }
  return response;
}

void DiagnoseSmoke303ImpulseResponse() {
  constexpr double sampleRate = 48'000.0;
  const auto response = RenderSmoke303Excitation();
  const auto report = [&](const char *name, const auto &response) {
    double worstPeriodicity = 0.0;
    for (const auto band : {std::array<double, 2>{80.0, 160.0},
                            std::array<double, 2>{160.0, 320.0},
                            std::array<double, 2>{320.0, 640.0},
                            std::array<double, 2>{640.0, 1'280.0},
                            std::array<double, 2>{1'280.0, 2'560.0}})
      worstPeriodicity = std::max(
          worstPeriodicity,
          LowFrequencyPeriodicity(response, sampleRate, band[0], band[1], 8.1));
    const double earlyTail = SegmentEnergy(response, 388'800, 432'000);
    const double lateTail = SegmentEnergy(response, 480'000, 528'000);
    const double gaussian =
        LateGaussianDeviation(response, 388'800, 480'000, 2'400);
    std::cout << "Smoke-303 topology " << name
              << ": worst periodicity=" << worstPeriodicity
              << " Gaussian deviation=" << gaussian
              << " late/early energy=" << lateTail / std::max(earlyTail, 1.e-20)
              << '\n';
  };
  report("fixed reference VFM", response);
}

} // namespace

int main() {
  TestDecayFilterRejectsInvalidT60();
  TestCompleteLatePathsFlushSubnormals();
  TestLateReverbAtEverySupportedRate();
  TestVelvetFeedbackMatrixIsParaunitaryAndDense();
  TestVelvetFractionalModulationPreservesTheStaticPath();
  TestLateReverbFlavourSelectionAndTransition();
  TestDiffusionControlsTemporalSpread();
  TestRoomScaleControlsTheCompleteVelvetSpan();
  TestRoomGeometryControlsLateArrivalTime();
  TestDampingProducesAudibleFrequencyDependentDecay();
  TestImpulseIsFiniteDenseAndStereo();
  TestDefaultBassTailAvoidsPeriodicEchoBuildUp();
  TestDecayControl();
  TestMidbandT60TracksDecayControl();
  TestMidbandT60RemainsIndependentOfRoomSize();
  TestMaximumRoomLongDecayAndModulationRemainStableWithoutShimmer();
  TestWindowedPitchShifterRaisesAnOctave();
  TestPitchShifterAdvancePreservesWarmState();
  TestPitchShifterHasNoGrainBoundarySpikes();
  TestPitchShifterRejectsAliasingInput();
  TestShimmerProjectionIsOrthonormal();
  TestFixedStereoFdnVectorsAreOrthonormal();
  TestShimmerReturnsOnlyThroughTheMainFeedbackLoop();
  TestShimmerAddsAnOctaveLayerWithoutMutingTheTail();
  TestShimmerCreatesAnOctaveInsideTheLateField();
  TestMaximumShimmerDoesNotGrowAfterExcitation();
  TestMaximumShimmerIsStableAtProductionSampleRate();
  TestShimmerLayerIsDiffuseAndNonPeriodic();
  TestSizeAutomationDoesNotDopplerTheLongTail();
  TestFixedStereoEncodingAndModulationRangeAffectResponse();
  TestAutomatedLateControlsPreserveStoredTail();
  DiagnoseSmoke303ImpulseResponse();
  std::cout << "Late reverb tests passed\n";
}
