import assert from "node:assert/strict";

import { calibrationParameterValues, calibrationPatch } from
  "../web/instrument_calibrations.mjs";

const descriptors = [
  ["impact_width", 1, .25, 4, "logarithmic"],
  ["bloom_energy_dependence", .35, 0, 1, "linear"],
  ["body_brightness", 0, -24, 24, "linear"],
  ["field_turbulence", .4, 0, 1, "linear"],
  ["body_decay_seconds_0", 3, .05, 20, "logarithmic"],
  ["body_decay_active_1", 0, 0, 1, "boolean"],
  ["body_decay_active_6", 0, 0, 1, "boolean"],
  ["resolved_frequency_0", 100, 40, 22000, "logarithmic"],
  ["resolved_level_0", 0, -24, 24, "linear"],
  ["model_level_db", -36, -60, 0, "linear"],
  ["impact_chirp_pitch", 1, .05, 4, "logarithmic"],
  ["body_decay_seconds_7", 1.2, .02, 20, "logarithmic"],
  ["body_low_cut", 40, 10, 1000, "logarithmic"],
  ["body_colour_gain", .5, -18, 18, "linear"],
].map(([key, defaultValue, minimum, maximum, scale], index) => ({
  key, defaultValue, minimum, maximum, scale, index,
}));

const defaults = calibrationParameterValues(
  {}, descriptors);
assert.deepEqual(defaults, descriptors.map(item => item.defaultValue));

assert.throws(() => calibrationParameterValues(
  { parameter_preset: "crash-start" }, descriptors), /Unknown/);

const gong = calibrationParameterValues(
  { parameter_preset: "gong-v1" }, descriptors);
assert.deepEqual(gong, [
  2.2, 1, -10, .72, 12, 0, 0, 128.9, -23, -35, .15, 1.7, 25, 2,
]);

const membraneDescriptors = [
  { key: "model_level_db", defaultValue: -10, index: 0 },
  { key: "fundamental_hz", defaultValue: 105, index: 1 },
  { key: "decay_seconds", defaultValue: 1.15, index: 2 },
  { key: "contact_direct_level", defaultValue: .245, index: 3 },
  { key: "contact_body_level", defaultValue: .7, index: 4 },
  { key: "fm_direct_level", defaultValue: .0144, index: 5 },
  { key: "fm_body_level", defaultValue: .081, index: 6 },
];
const acousticKick = calibrationParameterValues(
  { parameter_preset: "acoustic-kick" }, membraneDescriptors);
assert.deepEqual(acousticKick, [-12, 35, .25, 1.64, .205, .04, .05625]);
const acousticPatch = calibrationPatch(
  { parameter_preset: "acoustic-kick" }, membraneDescriptors,
  acousticKick, null,
);
assert.equal(acousticPatch.name, "Acoustic kick");
assert.equal(acousticPatch.recipe, "drum.membrane.v1");
assert.ok(acousticPatch.connections.every(connection =>
  !Object.hasOwn(connection, "gain")));
assert.equal(acousticPatch.nodes.find(
  node => node.id === "membrane-direct-mix")
  .parameters.contact_direct_level, 1.64);
assert.equal(acousticPatch.nodes.find(
  node => node.id === "membrane-body-mix").parameters.fm_body_level, .05625);

console.log("instrument calibration tests passed");
