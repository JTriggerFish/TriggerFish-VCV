import assert from "node:assert/strict";

import { calibrationParameterValues, calibrationPatch } from
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

const crash = calibrationParameterValues(
  { parameter_preset: "crash-start" }, descriptors);
assert.equal(crash[0], 1);
assert.equal(crash[2], 0);
assert.equal(crash[4], 8);
assert.equal(crash[6], -18);

const membraneDescriptors = [
  { key: "model_level_db", defaultValue: -10, index: 0 },
  { key: "fundamental_hz", defaultValue: 105, index: 1 },
  { key: "decay_seconds", defaultValue: 1.15, index: 2 },
];
const acousticKick = calibrationParameterValues(
  { parameter_preset: "acoustic-kick" }, membraneDescriptors);
assert.deepEqual(acousticKick, [-12, 35, .25]);
const acousticPatch = calibrationPatch(
  { parameter_preset: "acoustic-kick" }, membraneDescriptors,
  acousticKick, null,
);
assert.equal(acousticPatch.name, "Acoustic kick");
assert.equal(acousticPatch.recipe, "drum.membrane.v1");
assert.equal(acousticPatch.connections.find(
  item => item.id === "membrane-route-contact-direct").gain, 2);
assert.equal(acousticPatch.connections.find(
  item => item.id === "membrane-route-contact-body").gain, .25);
assert.equal(acousticPatch.connections.find(
  item => item.id === "membrane-route-fm-body").gain, .1125);

const ride = calibrationParameterValues(
  { parameter_preset: "ride-start" }, descriptors);
assert.equal(ride[0], .82);
assert.equal(ride[1], .28);
assert.equal(ride[2], .5);
assert.equal(ride[3], .52);
assert.equal(ride[4], 12);

const hihat = calibrationParameterValues(
  { parameter_preset: "hihat-start" }, descriptors);
assert.equal(hihat[0], .62);
assert.equal(hihat[1], .2);
assert.equal(hihat[2], 2.2);
assert.equal(hihat[3], .82);
assert.equal(hihat[4], 1.2);
assert.equal(hihat[5], 620);
assert.equal(hihat[6], -15);

const gong = calibrationParameterValues(
  { parameter_preset: "gong-start" }, descriptors);
assert.equal(gong[0], 3.1);
assert.equal(gong[1], .65);
assert.equal(gong[2], 1);
assert.equal(gong[5], 120);
assert.equal(gong[6], -8);

console.log("instrument calibration tests passed");
