#include "dsp_test_support.hpp"
#include "models/Tb303Voice.hpp"

#include <limits>

using dsp_test::Check;
using dsp_test::CheckNear;

int main() {
  for (const int rate : {44100, 48000, 96000, 192000}) {
    tfdsp::Tb303ResonanceControl control;
    control.SetSampleRate(rate);
    CheckNear(control.ProcessKnob(.8), .8, 0.0,
              "a loaded patch starts at its saved resonance");
    control.Reset();
    CheckNear(control.ProcessKnob(0.0), 0.0, 0.0,
              "reset initializes from the next parameter value");
    double value = 0.0;
    const int samples = static_cast<int>(std::lround(rate * .005));
    for (int i = 0; i < samples; ++i) {
      const double previous = value;
      value = control.ProcessKnob(1.0);
      Check(value > previous && value < 1.0,
            "a direct parameter jump is smoothed monotonically");
      for (int channel = 0; channel < 16; ++channel) {
        const double cv = (channel - 8) * .5;
        CheckNear(control.Modulate(value, -.75, cv), value - .075 * cv, 1.e-15,
                  "polyphonic CV remains immediate and attenuverted");
      }
    }
    CheckNear(value, 1.0 - std::exp(-samples / (rate * .005)), 1.e-12,
              "the smoothing time is consistent across sample rates");
    control.SetSampleRate(rate * 2);
    const double next = control.ProcessKnob(1.0);
    Check(next > value && next - value < .001,
          "changing sample rate preserves control continuity");
    control.Reset();
    CheckNear(control.ProcessKnob(.25), .25, 0.0,
              "reset does not retain the preceding automation target");
    CheckNear(
        control.Modulate(.25, 1.0, std::numeric_limits<double>::quiet_NaN()),
        .25, 0.0, "invalid CV does not contaminate knob state");
    Check(std::isfinite(
              control.ProcessKnob(std::numeric_limits<double>::infinity())),
          "invalid parameter updates remain finite");
  }
  return dsp_test::failures ? 1 : 0;
}
