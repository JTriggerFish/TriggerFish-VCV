// Pointer coordinates must include SVG letterboxing, not just its CSS bounds.
export function svgPosition(matrix, x, y) {
  if (!matrix) return null;
  const { a, b, c, d, e, f } = matrix;
  const determinant = a * d - b * c;
  if (!Number.isFinite(determinant) || Math.abs(determinant) < 1e-12) return null;
  return { x: (d * (x - e) - c * (y - f)) / determinant,
    y: (a * (y - f) - b * (x - e)) / determinant };
}

// A one-second knee gives long cymbal decays room without losing short settings.
// Stored values and interpolation remain log2(seconds); only the display changes.
export function decayPosition(seconds, minimum, maximum) {
  const low = Math.log1p(minimum), high = Math.log1p(maximum);
  return (Math.log1p(seconds) - low) / (high - low);
}

export function decaySeconds(position, minimum, maximum) {
  const low = Math.log1p(minimum), high = Math.log1p(maximum);
  return Math.expm1(low + position * (high - low));
}

export function decayDragDelta(previous, current, minimum, maximum, height, fine, anchor) {
  const low = 2 ** minimum, high = 2 ** maximum;
  const start = decayPosition(2 ** anchor, low, high);
  const delta = -(current - previous) / height * (fine ? .1 : 1);
  const end = Math.max(0, Math.min(1, start + delta));
  return Math.log2(decaySeconds(end, low, high)) - anchor;
}

export function shiftDecayPoints(points, delta, minimum, maximum) {
  const lower = minimum - Math.min(...points.map(point => point.y));
  const upper = maximum - Math.max(...points.map(point => point.y));
  const shift = Math.max(lower, Math.min(upper, delta));
  return points.map(point => ({ ...point, y: point.y + shift }));
}
