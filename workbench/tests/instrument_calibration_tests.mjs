import assert from "node:assert/strict";

import { calibrationParameterValues, calibrationPatch } from
  "../web/instrument_calibrations.mjs";

const descriptors = [
  ["impact_width", 1, .25, 4, "logarithmic"],
  ["bloom_energy_acceleration", .35, 0, 1, "linear"],
  ["body_brightness", 0, -72, 24, "linear"],
  ["body_excitation_centre", 1000, 40, 15000, "logarithmic"],
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
  2.2, .8, -56, 1100, .72, 12, 0, 0, 128.9, -12.12, 0, .15,
  1.1, 25, 2,
]);

const kickDescriptors = [
  {key:"model_level_db", defaultValue:-12, index:0},
  {key:"thump_pitch_hz", defaultValue:28, index:1},
  {key:"contact_level", defaultValue:.4, index:2},
];
const kickValues = calibrationParameterValues({parameter_preset:"kick"}, kickDescriptors);
const kickPatch = calibrationPatch({parameter_preset:"kick"}, kickDescriptors, kickValues, null);
assert.equal(kickPatch.recipe, "drum.kick.v1");
assert.equal(kickPatch.nodes.find(n=>n.id==="kick-contact").parameters.contact_level, .4);
assert.throws(()=>calibrationParameterValues({parameter_preset:"acoustic-kick"},kickDescriptors), /Unknown/);
console.log("instrument calibration tests passed");
