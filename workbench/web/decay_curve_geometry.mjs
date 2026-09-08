// Pointer coordinates must include SVG letterboxing, not just its CSS bounds.
export function svgPosition(matrix, x, y) {
  if (!matrix) return null;
  const { a, b, c, d, e, f } = matrix;
  const determinant = a * d - b * c;
  if (!Number.isFinite(determinant) || Math.abs(determinant) < 1e-12) return null;
  return { x: (d * (x - e) - c * (y - f)) / determinant,
    y: (a * (y - f) - b * (x - e)) / determinant };
}

export function decayDragDelta(previous, current, minimum, maximum, height, fine) {
  return -(current - previous) * (maximum - minimum) / height * (fine ? .1 : 1);
}

export function shiftDecayPoints(points, delta, minimum, maximum) {
  const lower = minimum - Math.min(...points.map(point => point.y));
  const upper = maximum - Math.max(...points.map(point => point.y));
  const shift = Math.max(lower, Math.min(upper, delta));
  return points.map(point => ({ ...point, y: point.y + shift }));
}
