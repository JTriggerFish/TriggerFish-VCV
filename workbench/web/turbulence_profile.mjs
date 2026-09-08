// Mirrors EvaluateTurbulence for editor geometry, not a second audio engine.
const clamp = (x, low, high) => Math.max(low, Math.min(high, Number.isFinite(x) ? x : 0));
export function turbulenceIntensity(frequency, level, slope, centre, local, relaxed) {
  const octaves = Math.log2(clamp(frequency, 1, 384000) / clamp(centre, 1, 384000));
  slope = clamp(slope, -1, 1);
  local = clamp(local, 0, 2);
  return relaxed ? clamp(level, 0, 4) * 2 ** (slope * octaves) * local
    : clamp(clamp(clamp(level, 0, 1) + slope * octaves, 0, 1) * local, 0, 1);
}
