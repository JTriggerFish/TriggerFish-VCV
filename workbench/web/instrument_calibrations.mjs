import {
  createAcousticKickPatch, membranePresetValues, patchWithMembraneValues,
} from "./membrane_patch.mjs";

const CrashStart = {
  model_level_db: -39,
  direct_gain: .05,
  impact_tone_noise: .96,
  impact_width: 1,
  field_gain: .00376872,
  body_brightness: 0,
  direct_low_cut: 100,
  direct_colour_frequency: 6500,
  direct_colour_gain: -6,
  direct_high_cut: 12000,
  dense_high_cut: 14000,
  body_decay_seconds_0: 8,
  body_decay_seconds_1: 7,
  body_decay_seconds_2: 6,
  body_decay_seconds_3: 4,
  body_decay_seconds_4: 1.5,
  body_decay_seconds_5: .8,
  body_decay_seconds_6: .4,
  body_decay_seconds_7: .02,
  body_decay_active_1: 1,
  body_decay_active_2: 1,
  body_decay_active_3: 1,
  body_decay_active_4: 1,
  body_decay_active_5: 1,
  body_decay_active_6: 1,
};

const CrashLevels = [
  -18, -15, -7, -5, -7, -7, -8, -5, 0, -4, -2, 0,
  -2, -3, -1, -3, -5, -6, -8, -12, -14, -15, -16, -18,
];
CrashLevels.forEach((level, index) => {
  CrashStart[`resolved_level_${index}`] = level;
});

const GongFrequencies = [
  120, 150, 180, 220, 300, 400, 520, 680, 850, 1100, 1400, 1750,
  2200, 2800, 3500, 4300, 5200, 6200, 7300, 8500, 9800, 11200,
  12800, 14500,
];
const GongLevels = [
  -8, 5, 6, 6, 5, 4, 3, 3, 2, 2, 1, 0,
  0, -1, -2, -3, -4, -5, -6, -8, -10, -12, -14, -16,
];
const GongStart = {
  model_level_db: -22.5,
  direct_gain: .2,
  impact_tone_noise: .3,
  impact_width: 3.1,
  bloom_level: 1.3,
  bloom_nonlinearity: .65,
  bloom_development: .85,
  bloom_diffusion: .8,
  field_gain: .018,
  body_brightness: 1,
  field_turbulence: .58,
  field_packet_spread: 4.5,
  field_phase_bandwidth: .5,
  field_exchange: .55,
  direct_low_cut: 50,
  direct_colour_frequency: 300,
  direct_colour_gain: 8,
  direct_high_cut: 6000,
  dense_low_cut: 25,
  dense_colour_frequency: 5000,
  dense_colour_gain: -1,
  dense_high_cut: 14000,
  body_decay_frequency_1: 150,
  body_decay_frequency_2: 500,
  body_decay_frequency_3: 1500,
  body_decay_frequency_4: 4000,
  body_decay_frequency_5: 8000,
  body_decay_frequency_6: 16000,
  body_decay_seconds_0: 3,
  body_decay_seconds_1: 6,
  body_decay_seconds_2: 8,
  body_decay_seconds_3: 10,
  body_decay_seconds_4: 9,
  body_decay_seconds_5: 7,
  body_decay_seconds_6: 3,
  body_decay_seconds_7: .2,
  body_decay_active_1: 1,
  body_decay_active_2: 1,
  body_decay_active_3: 1,
  body_decay_active_4: 1,
  body_decay_active_5: 1,
  body_decay_active_6: 1,
};
GongFrequencies.forEach((frequency, index) => {
  GongStart[`resolved_frequency_${index}`] = frequency;
  GongStart[`resolved_level_${index}`] = GongLevels[index];
});

const RideStart = {
  model_level_db: -18.3,
  impact_tone_noise: .72,
  impact_width: .82,
  bloom_level: .42,
  bloom_nonlinearity: .28,
  bloom_development: .3,
  bloom_diffusion: .38,
  body_brightness: .5,
  field_turbulence: .52,
  field_packet_spread: 2.8,
  field_phase_bandwidth: .24,
  field_exchange: .55,
  field_gain: .0376872,
  direct_gain: .5,
  direct_colour_frequency: 3000,
  direct_colour_gain: -12,
  direct_high_cut: 1200,
  body_decay_seconds_0: 12,
  body_decay_seconds_1: 9,
  body_decay_seconds_2: 7,
  body_decay_seconds_3: 1.5,
  body_decay_seconds_4: 1,
  body_decay_seconds_5: .5,
  body_decay_seconds_6: .3,
  body_decay_seconds_7: .25,
  body_decay_active_1: 1,
  body_decay_active_2: 1,
  body_decay_active_3: 1,
  body_decay_active_4: 1,
  body_decay_active_5: 1,
  body_decay_active_6: 1,
};

const HiHatFrequencies = [
  620, 760, 930, 1120, 1350, 1630, 1960, 2350, 2800, 3300, 3850, 4450,
  5100, 5800, 6550, 7350, 8200, 9100, 10000, 10900, 11800, 12700,
  13800, 15000,
];
const HiHatLevels = [
  -15, -12, -10, -8, -7, -6, -5, -4, -3, -2, -1, 0,
  1, 2, 3, 3, 4, 4, 3, 2, 1, 0, -2, -5,
];
const HiHatStart = {
  model_level_db: -16.5,
  direct_gain: .08,
  impact_tone_noise: .94,
  impact_width: .62,
  bloom_level: .08,
  bloom_nonlinearity: .2,
  bloom_development: .14,
  bloom_diffusion: .3,
  field_gain: .12,
  body_brightness: 2.2,
  field_turbulence: .82,
  field_packet_spread: 5.5,
  field_phase_bandwidth: .8,
  field_exchange: .72,
  direct_low_cut: 420,
  direct_colour_frequency: 7200,
  direct_colour_gain: 4,
  direct_high_cut: 19000,
  dense_low_cut: 320,
  dense_colour_frequency: 8000,
  dense_colour_gain: 3,
  dense_high_cut: 19000,
  body_decay_frequency_1: 500,
  body_decay_frequency_2: 1500,
  body_decay_frequency_3: 4000,
  body_decay_frequency_4: 8000,
  body_decay_frequency_5: 12000,
  body_decay_frequency_6: 17000,
  body_decay_seconds_0: 1.2,
  body_decay_seconds_1: 1.7,
  body_decay_seconds_2: 2.1,
  body_decay_seconds_3: 2.2,
  body_decay_seconds_4: 1.7,
  body_decay_seconds_5: 1.15,
  body_decay_seconds_6: .65,
  body_decay_seconds_7: .08,
  body_decay_active_1: 1,
  body_decay_active_2: 1,
  body_decay_active_3: 1,
  body_decay_active_4: 1,
  body_decay_active_5: 1,
  body_decay_active_6: 1,
};
HiHatFrequencies.forEach((frequency, index) => {
  HiHatStart[`resolved_frequency_${index}`] = frequency;
  HiHatStart[`resolved_level_${index}`] = HiHatLevels[index];
});

function applyOverrides(values, descriptors, overrides) {
  for (const descriptor of descriptors) {
    if (Object.hasOwn(overrides, descriptor.key))
      values[descriptor.index] = overrides[descriptor.key];
  }
}

export function calibrationParameterValues(calibration, descriptors) {
  if (calibration.parameter_preset === "acoustic-kick")
    return membranePresetValues("acousticKick", descriptors);
  const values = descriptors.map(item => item.defaultValue);
  if (calibration.parameter_preset === "crash-start") {
    applyOverrides(values, descriptors, CrashStart);
  } else if (calibration.parameter_preset === "hihat-start") {
    applyOverrides(values, descriptors, HiHatStart);
  } else if (calibration.parameter_preset === "gong-start") {
    applyOverrides(values, descriptors, GongStart);
  } else if (calibration.parameter_preset === "ride-start") {
    applyOverrides(values, descriptors, RideStart);
  }
  return values;
}

export function calibrationPatch(
  calibration, descriptors, values, fallbackPatch,
) {
  if (calibration.parameter_preset !== "acoustic-kick") return fallbackPatch;
  return patchWithMembraneValues(
    createAcousticKickPatch(descriptors), descriptors, values,
  );
}
