const clamp = (value, minimum, maximum) =>
  Math.max(minimum, Math.min(maximum, value));

const gongPreset = {
  impact_tone_noise: .35, impact_width: 1.8,
  bloom_level: 1.15, bloom_nonlinearity: .65, bloom_development: .8,
  body_tone_wash: .48, body_brightness: -2.5,
  dense_minimum_frequency: 120, dense_maximum_frequency: 14000,
  dense_mode_density: 1,
  dense_spacing_jitter: .72, dense_gain_spread: 3,
  dense_decay_spread: .25,
  turbulence_amount: .15, turbulence_persistence: 1.3,
  turbulence_frequency_0: 120, turbulence_frequency_1: 800,
  turbulence_frequency_2: 5000,
  body_decay_frequency_1: 240, body_decay_frequency_2: 700,
  body_decay_frequency_3: 2400, body_decay_frequency_4: 8000,
  body_decay_seconds_0: 12, body_decay_seconds_1: 10,
  body_decay_seconds_2: 7, body_decay_seconds_3: 4,
  body_decay_seconds_4: 3, body_decay_seconds_7: 2.5,
  body_decay_active_1: 1, body_decay_active_2: 1,
  body_decay_active_3: 1, body_decay_active_4: 1,
};

const gongLevels = [5, 6, 8, 9, 8, 7, 5, 3, 1, -1, -3, -5];
const gongWashLevels = [7, 8, 6, 3, 0, -2, -4, -6];

const chinaPreset = {
  impact_tone_noise: .62, impact_width: .65,
  bloom_level: .85, bloom_nonlinearity: .5, bloom_development: .3,
  body_tone_wash: .7, body_brightness: 1.5,
  dense_minimum_frequency: 350, dense_maximum_frequency: 22000,
  dense_mode_density: .65,
  dense_spacing_jitter: .92, dense_gain_spread: 5,
  dense_decay_spread: .1,
  turbulence_amount: .12, turbulence_persistence: .65,
  turbulence_frequency_0: 900, turbulence_frequency_1: 5000,
  turbulence_frequency_2: 16000,
  body_decay_frequency_1: 1100, body_decay_frequency_2: 3000,
  body_decay_frequency_3: 8000, body_decay_frequency_4: 17000,
  body_decay_seconds_0: 2.8, body_decay_seconds_1: 2.2,
  body_decay_seconds_2: 1.5, body_decay_seconds_3: .8,
  body_decay_seconds_4: .5, body_decay_seconds_7: .35,
  body_decay_active_1: 1, body_decay_active_2: 1,
  body_decay_active_3: 1, body_decay_active_4: 1,
};

const chinaLevels = [-7, -5, -3, 0, 4, 7, 6, 5, 3, 2, 0, -2];
const chinaWashLevels = [-7, -4, 0, 5, 6, 3, 0, -2];

function presetWithSpectralControls(
  preset, defaults, scale, resolvedLevels, washLevels,
) {
  const result = { ...preset };
  for (let index = 0; index < 12; ++index) {
    result[`resolved_frequency_${index}`] =
      defaults[`resolved_frequency_${index}`] * scale;
    result[`resolved_level_${index}`] = resolvedLevels[index];
  }
  for (let index = 0; index < 8; ++index) {
    result[`dense_wash_frequency_${index}`] =
      defaults[`dense_wash_frequency_${index}`] * scale;
    result[`dense_wash_level_${index}`] = washLevels[index];
  }
  return result;
}

export function expandedSizeMeta(descriptors, position) {
  const amount = clamp(Number(position), 0, 1);
  const defaults = Object.fromEntries(
    descriptors.map(item => [item.key, item.defaultValue]),
  );
  const gong = presetWithSpectralControls(
    gongPreset, defaults, .42, gongLevels, gongWashLevels,
  );
  const china = presetWithSpectralControls(
    chinaPreset, defaults, 1.45, chinaLevels, chinaWashLevels,
  );
  const target = amount < .5 ? gong : china;
  const blend = amount < .5 ? amount * 2 : (amount - .5) * 2;
  return Object.entries(target).flatMap(([key, endpoint]) => {
    const descriptor = descriptors.find(item => item.key === key);
    if (!descriptor) return [];
    const neutral = defaults[key];
    const left = amount < .5 ? endpoint : neutral;
    const right = amount < .5 ? neutral : endpoint;
    const value = descriptor.scale === "logarithmic"
      ? left * (right / left) ** blend
      : left + (right - left) * blend;
    return [{ descriptor, value: clamp(
      value, descriptor.minimum, descriptor.maximum,
    ) }];
  });
}
