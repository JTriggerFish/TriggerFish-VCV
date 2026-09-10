// Shared by the actual slider and the timing popover's read-only preview.
const SlowestBloomRate = .01;
export const bloomRateNormalized = (descriptor, value) => value <= 0 ? 0 :
  Math.max(0, Math.min(1, .02 + .98 * Math.log(value / SlowestBloomRate) /
    Math.log(descriptor.maximum / SlowestBloomRate)));
export const bloomRateDenormalized = (descriptor, position) => position < .01
  ? 0 : SlowestBloomRate * (descriptor.maximum / SlowestBloomRate) **
    ((position - .02) / .98);
