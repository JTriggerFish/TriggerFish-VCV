// UI preview of modal_packet_allocator.hpp; no audio synthesis here.
export function packetAllocation(points, density, width, pairedCentres = false) {
  const active = points.filter(p => p.active).length;
  const centresPerHandle = pairedCentres ? 2 : 1;
  const centres = active * centresPerHandle;
  const pairs = points.map(() => 0);
  const weights = points.map(p => p.active && width(p) > 0
    ? Math.sqrt(width(p)) * Math.max(0, Math.min(4, p.allocation ?? 1)) : 0);
  const total = weights.reduce((a, b) => a + b, 0);
  const budget = Math.round(Math.floor((512 - centres) / 2) * density);
  if (!(total > 0)) return {pairs, count: centres, centresPerHandle};
  const fractions = weights.map((weight, i) => {
    const exact = budget * weight / total;
    pairs[i] = Math.floor(exact);
    return exact - pairs[i];
  });
  for (let left = budget - pairs.reduce((a, b) => a + b, 0); left > 0; --left) {
    const best = fractions.indexOf(Math.max(...fractions));
    pairs[best]++; fractions[best] = 0;
  }
  return {pairs, count: centres + 2 * budget, centresPerHandle};
}
