#pragma once

#include "enveloped_noise_burst.hpp"
#include "finite_force_pulse.hpp"
#include "micro_contact_burst.hpp"
#include "tonal_contact_chirp.hpp"
#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <cmath>

namespace tfdsp::percussion {

// Derived port projection supplied by a visible performance macro such as
// Implement. It is not serialized as independent patch state.
struct ContactPortProjection {
  float pulseDirect{0.f};
  float pulseBody{1.f};
  float chirpDirect{1.f};
  float chirpBody{.25f};
  float noiseDirect{.25f};
  float noiseBody{.5f};
  float microDirect{.5f};
  float microBody{.75f};
};

struct ContactExciterParameters {
  float pulseDurationSeconds{.001f};
  float pulseAmplitude{1.f};
  TonalContactChirpParameters chirp{};
  EnvelopedNoiseBurstParameters noise{};
  MicroContactBurstParameters microContacts{};
  ContactPortProjection projection{};
};

struct ContactExciterSample {
  float directRadiation{};
  float bodyDrive{};
};

// Composes independent contact primitives into explicit direct-radiation and
// body-drive ports. Implement and hardness macros select parameters upstream.
class ContactExciter {
public:
  void Prepare(const float sampleRate) {
    pulse_.Prepare(sampleRate);
    chirp_.Prepare(sampleRate);
    noise_.Prepare(sampleRate);
    microContacts_.Prepare(sampleRate);
    Reset();
  }

  void Reset() noexcept {
    pulse_.Reset();
    chirp_.Reset();
    noise_.Reset();
    microContacts_.Reset();
  }

  void Trigger(const ContactExciterParameters &parameters) noexcept {
    projection_ = parameters.projection;
    SanitizeProjection();
    pulse_.Trigger(parameters.pulseDurationSeconds, parameters.pulseAmplitude);
    chirp_.Trigger(parameters.chirp);
    noise_.Trigger(parameters.noise);
    microContacts_.Trigger(parameters.microContacts);
  }

  ContactExciterSample Process() noexcept {
    const float pulse = pulse_.Process();
    const float chirp = chirp_.Process();
    const float noise = noise_.Process();
    const float micro = microContacts_.Process();
    ContactExciterSample output;
    output.directRadiation = projection_.pulseDirect * pulse +
        projection_.chirpDirect * chirp + projection_.noiseDirect * noise +
        projection_.microDirect * micro;
    output.bodyDrive = projection_.pulseBody * pulse +
        projection_.chirpBody * chirp + projection_.noiseBody * noise +
        projection_.microBody * micro;
    if (!std::isfinite(output.directRadiation))
      output.directRadiation = 0.f;
    if (!std::isfinite(output.bodyDrive))
      output.bodyDrive = 0.f;
    return output;
  }

  bool Active() const noexcept {
    return pulse_.Active() || chirp_.Active() || noise_.Active() ||
           microContacts_.Active();
  }

private:
  static float BoundedGain(const float gain) noexcept {
    return std::clamp(tfdsp::FiniteNormalOrZero(gain), -16.f, 16.f);
  }

  void SanitizeProjection() noexcept {
    projection_.pulseDirect = BoundedGain(projection_.pulseDirect);
    projection_.pulseBody = BoundedGain(projection_.pulseBody);
    projection_.chirpDirect = BoundedGain(projection_.chirpDirect);
    projection_.chirpBody = BoundedGain(projection_.chirpBody);
    projection_.noiseDirect = BoundedGain(projection_.noiseDirect);
    projection_.noiseBody = BoundedGain(projection_.noiseBody);
    projection_.microDirect = BoundedGain(projection_.microDirect);
    projection_.microBody = BoundedGain(projection_.microBody);
  }

  FiniteForcePulse pulse_{};
  TonalContactChirp chirp_{};
  EnvelopedNoiseBurst noise_{};
  MicroContactBurst microContacts_{};
  ContactPortProjection projection_{};
};

} // namespace tfdsp::percussion
