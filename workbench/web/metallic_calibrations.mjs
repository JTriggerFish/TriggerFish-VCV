const GongFrequencies = [
  128.9, 245, 304.7, 375, 421.9, 550.8, 621.1, 726.6,
  878.9, 1380, 1699.2, 1980.5, 2900, 4200, 6000, 8800, 12000,
];

const GongLevels = [
  -12.12, -22.62, -16.47, -8.8, -16.09, -7.95, -15.81,
  -13.31, -16.1, -5.5, -5.5, -7.3, -5.5, -4, -4, 6, 0,
];

const GongTurbulence = GongFrequencies.map(frequency =>
  frequency < 1000 ? .12 : frequency < 1500 ? .3 :
    frequency < 1800 ? .4 : frequency < 2500 ? .5 : 1.2);

const GongV1 = {
  // Gain-unit reset only: this remains an unvalidated timbral starting point.
  model_level_db: 0,
  impact_tone_noise: .2,
  impact_width: 2.2,
  impact_chirp_pitch: .15,
  bloom_rate: .5,
  bloom_energy_acceleration: .8,
  bloom_phase_diffusion: .8,
  body_brightness: -56,
  body_excitation_centre: 1100,
  body_excitation: 1,
  field_turbulence: .72,
  field_turbulence_slope: .6,
  field_turbulence_centre: 1200,
  field_packet_spread: 5.5,
  field_satellite_density: .65,
  field_phase_bandwidth: .8,
  field_exchange: .5,
  field_gain: 4,
  direct_gain: .008,
  direct_high_cut: 3000,
  body_low_cut: 25,
  body_colour_frequency: 8500,
  body_colour_gain: 2,
  body_high_cut: 12000,
  body_decay_seconds_0: 12,
  body_decay_seconds_7: 1.1,
  velocity_brightness: 5,
};

for (let index = 1; index <= 6; ++index)
  GongV1[`body_decay_active_${index}`] = 0;

for (let index = 0; index < 32; ++index) {
  GongV1[`resolved_frequency_${index}`] = 15000;
  GongV1[`resolved_level_${index}`] = -72;
  GongV1[`resolved_turbulence_${index}`] = 1;
}

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
