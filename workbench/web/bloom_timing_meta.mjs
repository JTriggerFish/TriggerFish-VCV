// Design-time gesture only. Saves contain the expanded ordinary parameters.
export const BloomTimingKeys = Object.freeze([
  "bloom_rate", "body_brightness", "body_excitation_centre",
]);

export function bloomTimingValues(baseline, position, descriptors) {
  if (!Number.isFinite(position) || position < -1 || position > 1)
    throw Error("Bloom timing must be between -1 and 1");
  const targets = {
    bloom_rate: baseline.bloom_rate * 2 ** (-position),
    body_brightness: baseline.body_brightness - 6 * position,
    body_excitation_centre: baseline.body_excitation_centre * 2 ** (-.25 * position),
  };
  const values = {}, limited = [];
  for (const key of BloomTimingKeys) {
    const descriptor = descriptors.find(item => item.key === key);
    if (!descriptor || !Number.isFinite(baseline[key]) ||
        baseline[key] < descriptor.minimum || baseline[key] > descriptor.maximum)
      throw Error(`Invalid bloom timing baseline: ${key}`);
    values[key] = Math.max(descriptor.minimum, Math.min(descriptor.maximum, targets[key]));
    if (values[key] !== targets[key]) limited.push(descriptor.name ?? key);
  }
  return {values, limited};
}
