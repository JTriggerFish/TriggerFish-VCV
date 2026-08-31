#pragma once

#include "enveloped_noise_burst.hpp"
#include "finite_force_pulse.hpp"
#include "micro_contact_burst.hpp"
#include "tonal_contact_chirp.hpp"

#include <cmath>

namespace tfdsp::percussion {

struct ContactRoutingGains {
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
  ContactRoutingGains routing{};
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
    routing_ = parameters.routing;
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
    output.directRadiation = routing_.pulseDirect * pulse +
        routing_.chirpDirect * chirp + routing_.noiseDirect * noise +
        routing_.microDirect * micro;
    output.bodyDrive = routing_.pulseBody * pulse +
        routing_.chirpBody * chirp + routing_.noiseBody * noise +
        routing_.microBody * micro;
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
  ContactRoutingGains routing_{};
};

} // namespace tfdsp::percussion
