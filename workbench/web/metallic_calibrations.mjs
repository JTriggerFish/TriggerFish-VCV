const GongFrequencies = [
  128.9, 245, 304.7, 375, 421.9, 550.8, 621.1, 668,
  726.6, 878.9, 1100, 1380, 1699.2, 1980.5, 2400, 2900,
  3500, 4200, 5000, 6000, 7200, 8800, 11000, 15000,
];

const GongLevels = [
  -23, -17, -10, 0, -14, 4, -12.8, -17.3,
  -13.7, -14.3, -10, -8, -6, -6, -5, -4,
  -3, -2, -1, 0, 0, -2, -6, -12,
];

const GongTurbulence = GongFrequencies.map(frequency =>
  frequency < 500 ? .12 : frequency < 3500 ? .3 : 1.2);

const GongV1 = {
  model_level_db: -35,
  impact_tone_noise: .2,
  impact_width: 2.2,
  impact_chirp_pitch: .15,
  bloom_rate: 5,
  bloom_energy_dependence: 1,
  bloom_phase_diffusion: .8,
  body_brightness: -10,
  body_tilt_centre: 550,
  body_excitation: 1,
  field_turbulence: .72,
  field_turbulence_slope: .22,
  field_turbulence_centre: 1200,
  field_packet_spread: 5.5,
  field_phase_bandwidth: .8,
  field_exchange: .5,
  field_gain: .08,
  direct_gain: .002,
  direct_high_cut: 3000,
  body_low_cut: 25,
  body_colour_frequency: 8500,
  body_colour_gain: 2,
  body_high_cut: 12000,
  body_decay_seconds_0: 12,
  body_decay_seconds_7: 1.7,
  velocity_brightness: 5,
};

for (let index = 1; index <= 6; ++index)
  GongV1[`body_decay_active_${index}`] = 0;

GongFrequencies.forEach((frequency, index) => {
  GongV1[`resolved_frequency_${index}`] = frequency;
  GongV1[`resolved_level_${index}`] = GongLevels[index];
  GongV1[`resolved_turbulence_${index}`] = GongTurbulence[index];
});

export function metallicCalibrationValues(name, descriptors) {
  if (name !== "gong-v1") return null;
  return descriptors.map(descriptor => Object.hasOwn(GongV1, descriptor.key)
    ? GongV1[descriptor.key] : descriptor.defaultValue);
}

export function metallicCalibrationOverrides(name) {
  return name === "gong-v1" ? { ...GongV1 } : null;
}
