import assert from "node:assert/strict";

import { calibrationParameterValues, calibrationPatch, calibrationEvent, recipeStartingValues, recipeStartingEvent } from
  "../web/instrument_calibrations.mjs";
import KickCalibration from "../web/kick_calibration.fit.json" with { type: "json" };
import {referenceCalibration, checkedCalibrationValues} from "../web/reference_calibration_library.mjs";

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
  ["output_low_cut", 40, 10, 1000, "logarithmic"],
  ["output_colour_gain", .5, -18, 18, "linear"],
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
  2.2, .8, -56, 1100, .72*(1000/1200)**.6, 12, 0, 0, 128.9, -12.12, 0, .15,
  1.1, 25, 2,
]);

const kickParameters = Object.assign({}, ...KickCalibration.instrument.nodes.map(n=>n.parameters));
const kickDescriptors = Object.keys(kickParameters).map((key,index)=>({
  key,index,defaultValue:0,minimum:-Infinity,maximum:Infinity,
}));
const kickValues = calibrationParameterValues({parameter_preset:"kick"}, kickDescriptors);
const kickPatch = calibrationPatch({parameter_preset:"kick"}, kickDescriptors, kickValues, null);
assert.equal(kickPatch.recipe, "drum.kick.v1");
assert.equal(kickValues.length, 65);
assert.deepEqual(kickValues, Object.values(kickParameters));
assert.deepEqual(recipeStartingValues("drum.kick.v1",kickDescriptors),kickValues);
assert.deepEqual(recipeStartingEvent("drum.kick.v1"),KickCalibration.controls.event);
assert.deepEqual(kickPatch, KickCalibration.instrument);
assert.throws(()=>calibrationParameterValues({parameter_preset:"kick"},kickDescriptors.slice(1)), /surface/);
const invalidRange=kickDescriptors.map(d=>({...d,minimum:100000}));
assert.throws(()=>calibrationParameterValues({parameter_preset:"kick"},invalidRange), /Invalid/);
assert.throws(()=>calibrationParameterValues({parameter_preset:"acoustic-kick"},kickDescriptors), /Unknown/);

const snareFit = referenceCalibration("snare-standard");
const snareParameters = Object.assign({}, ...snareFit.instrument.nodes.map(n=>n.parameters));
const snareDescriptors = Object.keys(snareParameters).map((key,index)=>({
  key,index,defaultValue:0,minimum:-Infinity,maximum:Infinity,
}));
const snareTarget = {id:"snare-standard",recipe:"drum.snare.v1"};
const snareValues = calibrationParameterValues(snareTarget,snareDescriptors);
assert.deepEqual(snareValues,Object.values(snareParameters));
assert.deepEqual(recipeStartingValues("drum.snare.v1",snareDescriptors),snareValues);
assert.deepEqual(recipeStartingEvent("drum.snare.v1"),snareFit.controls.event);
assert.deepEqual(calibrationPatch(snareTarget,snareDescriptors,snareValues,null),snareFit.instrument);
assert.throws(()=>calibrationParameterValues({...snareTarget,recipe:"drum.kick.v1"},snareDescriptors),/recipe/);
assert.throws(()=>checkedCalibrationValues(snareFit,snareDescriptors.slice(1)),/surface/);
assert.throws(()=>checkedCalibrationValues(snareFit,snareDescriptors.map(d=>({...d,minimum:1e9}))),/Invalid/);
const duplicateFit = structuredClone(snareFit);
duplicateFit.instrument.nodes.push(structuredClone(duplicateFit.instrument.nodes[0]));
assert.throws(()=>checkedCalibrationValues(duplicateFit,snareDescriptors),/Duplicate/);
assert.equal(referenceCalibration("unknown"),null);

for (const id of ["crash-standard", "ride-standard", "gong-standard", "hihat-standard"]) {
  const fit = referenceCalibration(id);
  const parameters = Object.assign({}, ...fit.instrument.nodes.map(n=>n.parameters));
  const surface = Object.keys(parameters).map((key,index)=>({
    key,index,defaultValue:0,minimum:-Infinity,maximum:Infinity,
  }));
  const target = {id,recipe:fit.instrument.recipe};
  const values = calibrationParameterValues(target,surface);
  assert.deepEqual(values,Object.values(parameters));
  assert.deepEqual(calibrationPatch(target,surface,values,null),fit.instrument);
  // Reviewed shared curves may have interior points; this is not a DSP limit.
  const reviewedKnots = id === "crash-standard" ? []
    : id === "ride-standard" ? [700] : [];
  for (let knot=1;knot<=6;++knot) {
    const active = knot <= reviewedKnots.length;
    assert.equal(parameters[`body_decay_active_${knot}`],Number(active),
      `${id}: shared damping curve must match the reviewed preset`);
    if (active) assert.equal(parameters[`body_decay_frequency_${knot}`],reviewedKnots[knot-1]);
  }
}
const rideEvent = calibrationEvent({id:"ride-standard"});
assert.equal(rideEvent.strength, 0.5506513756876121);
assert.equal(rideEvent.location, 0.4904862153154393);
assert.equal(rideEvent.seed, 1944);
rideEvent.strength = 0;
assert.notEqual(calibrationEvent({id:"ride-standard"}).strength, 0,
  "loading an event must not mutate the saved calibration");
assert.equal(calibrationEvent({}), null);
console.log("instrument calibration tests passed");
