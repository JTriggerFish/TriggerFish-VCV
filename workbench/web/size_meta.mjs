import { metallicCalibrationOverrides } from "./metallic_calibrations.mjs";

const clamp = (value, minimum, maximum) =>
  Math.max(minimum, Math.min(maximum, value));

const chinaPreset = {
  impact_tone_noise: .62, impact_width: .65,
  bloom_rate: 5, bloom_energy_acceleration: .55,
  bloom_phase_diffusion: .8, body_brightness: 1.5,
  body_excitation_centre: 3000,
  field_turbulence_slope: .12,
  field_turbulence_centre: 4000,
  body_decay_seconds_0: 2.8, body_decay_seconds_7: .35,
  body_decay_active_1: 0, body_decay_active_2: 0,
  body_decay_active_3: 0, body_decay_active_4: 0,
  body_decay_active_5: 0, body_decay_active_6: 0,
};

const chinaLevels = [-7, -5, -3, 0, 4, 7, 6, 5, 3, 2, 0, -2];
function presetWithSpectralControls(
  preset, defaults, scale, resolvedLevels,
) {
  const result = { ...preset };
  const indices = Object.keys(defaults).flatMap(key => {
    const match = /^resolved_frequency_(\d+)$/.exec(key);
    return match ? [Number(match[1])] : [];
  }).sort((left, right) => left - right);
  indices.forEach((index, position) => {
    const curvePosition = indices.length > 1
      ? position * (resolvedLevels.length - 1) / (indices.length - 1) : 0;
    const left = Math.floor(curvePosition);
    const right = Math.min(resolvedLevels.length - 1, left + 1);
    const amount = curvePosition - left;
    result[`resolved_frequency_${index}`] =
      defaults[`resolved_frequency_${index}`] * scale;
    result[`resolved_level_${index}`] = resolvedLevels[left] + amount *
      (resolvedLevels[right] - resolvedLevels[left]);
  });
  return result;
}

export function expandedSizeMeta(descriptors, position) {
  const amount = clamp(Number(position), 0, 1);
  const defaults = Object.fromEntries(
    descriptors.map(item => [item.key, item.defaultValue]),
  );
  const gong = metallicCalibrationOverrides("gong-v1");
  const china = presetWithSpectralControls(
    chinaPreset, defaults, 1.45, chinaLevels,
  );
  const target = amount < .5 ? gong : china;
  const blend = amount < .5 ? amount * 2 : (amount - .5) * 2;
  const keys = new Set([...Object.keys(gong), ...Object.keys(china)]);
  return [...keys].flatMap(key => {
    const descriptor = descriptors.find(item => item.key === key);
    if (!descriptor) return [];
    const neutral = defaults[key];
    const endpoint = Object.hasOwn(target, key) ? target[key] : neutral;
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
