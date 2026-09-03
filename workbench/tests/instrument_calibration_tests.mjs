import assert from "node:assert/strict";

import { calibrationParameterValues } from
  "../web/instrument_calibrations.mjs";

const descriptors = [
  ["impact_width", 1, .25, 4, "logarithmic"],
  ["bloom_nonlinearity", .35, 0, 1, "linear"],
  ["body_brightness", 0, -4, 4, "linear"],
  ["field_turbulence", .4, 0, 1, "linear"],
  ["body_decay_seconds_0", 3, .05, 20, "logarithmic"],
  ["resolved_frequency_0", 100, 40, 22000, "logarithmic"],
  ["resolved_level_0", 0, -24, 24, "linear"],
].map(([key, defaultValue, minimum, maximum, scale], index) => ({
  key, defaultValue, minimum, maximum, scale, index,
}));

const defaults = calibrationParameterValues(
  { parameter_preset: "snare-default" }, descriptors);
assert.deepEqual(defaults, descriptors.map(item => item.defaultValue));

const ride = calibrationParameterValues(
  { parameter_preset: "ride-start" }, descriptors);
assert.equal(ride[0], .82);
assert.equal(ride[1], .28);
assert.equal(ride[2], -1.4);
assert.equal(ride[3], .52);
assert.equal(ride[4], 12);

const gong = calibrationParameterValues(
  { parameter_preset: "gong-start" }, descriptors);
assert.equal(gong[0], 1.8);
assert.equal(gong[1], .65);
assert.equal(gong[2], -2.5);
assert.equal(gong[5], 42);
assert.equal(gong[6], 5);

console.log("instrument calibration tests passed");
