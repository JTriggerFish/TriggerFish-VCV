#pragma once

#include "enveloped_noise_burst.hpp"
#include "finite_force_pulse.hpp"
#include "micro_contact_burst.hpp"
#include "tonal_contact_chirp.hpp"
#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <cmath>

namespace tfdsp::percussion {

struct ContactExciterParameters {
  float pulseDurationSeconds{.001f};
  float pulseAmplitude{1.f};
  TonalContactChirpParameters chirp{};
  EnvelopedNoiseBurstParameters noise{};
  MicroContactBurstParameters microContacts{};
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
    pulse_.Trigger(parameters.pulseDurationSeconds, parameters.pulseAmplitude);
    chirp_.Trigger(parameters.chirp);
    noise_.Trigger(parameters.noise);
    microContacts_.Trigger(parameters.microContacts);
    // Prepare unit conversions once per event, never from a playing tail.
    pulseBodyScale_ = pulse_.BodyImpulseScale();
    chirpBodyScale_ = chirp_.BodyImpulseScale();
    noiseBodyScale_ = noise_.BodyImpulseScale();
    microBodyScale_ = microContacts_.BodyImpulseScale();
  }

  ContactExciterSample Process() noexcept {
    const float pulse = pulse_.Process();
    const float chirp = chirp_.Process();
    const float noise = noise_.Process();
    const float micro = microContacts_.Process();
    ContactExciterSample output;
    output.directRadiation = pulse + chirp + noise + micro;
    output.bodyDrive = pulseBodyScale_ * pulse + chirpBodyScale_ * chirp +
        noiseBodyScale_ * noise + microBodyScale_ * micro;
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
  FiniteForcePulse pulse_{};
  TonalContactChirp chirp_{};
  EnvelopedNoiseBurst noise_{};
  MicroContactBurst microContacts_{};
  float pulseBodyScale_{};
  float chirpBodyScale_{};
  float noiseBodyScale_{};
  float microBodyScale_{};
};

} // namespace tfdsp::percussion
