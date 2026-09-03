import { expandedSizeMeta } from "./size_meta.mjs";

const RideStart = {
  impact_tone_noise: .72,
  impact_width: .82,
  bloom_level: .42,
  bloom_nonlinearity: .28,
  bloom_development: .3,
  bloom_diffusion: .38,
  body_brightness: -1.4,
  field_turbulence: .52,
  field_packet_spread: 2.8,
  field_phase_bandwidth: .24,
  field_exchange: .55,
  direct_gain: .34,
  body_decay_seconds_0: 12,
  body_decay_seconds_1: 9,
  body_decay_seconds_2: 7,
  body_decay_seconds_3: 3.5,
  body_decay_seconds_4: 2.2,
  body_decay_seconds_5: 1.2,
  body_decay_seconds_6: .7,
  body_decay_seconds_7: .25,
  body_decay_active_1: 1,
  body_decay_active_2: 1,
  body_decay_active_3: 1,
  body_decay_active_4: 1,
  body_decay_active_5: 1,
  body_decay_active_6: 1,
};

function applyOverrides(values, descriptors, overrides) {
  for (const descriptor of descriptors) {
    if (Object.hasOwn(overrides, descriptor.key))
      values[descriptor.index] = overrides[descriptor.key];
  }
}

export function calibrationParameterValues(calibration, descriptors) {
  const values = descriptors.map(item => item.defaultValue);
  if (calibration.parameter_preset === "gong-start") {
    for (const { descriptor, value } of expandedSizeMeta(descriptors, 0))
      values[descriptor.index] = value;
  } else if (calibration.parameter_preset === "ride-start") {
    applyOverrides(values, descriptors, RideStart);
  }
  return values;
}
