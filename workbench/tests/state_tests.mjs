import assert from "node:assert/strict";

import {
  fitMacroValues, readFit, snapshotState, validateFit,
} from "../web/state.mjs";
import { createKickPatch } from "../web/kick_patch.mjs";
import { createSnarePatch } from "../web/snare_patch.mjs";

const descriptors = [
  {
    index: 0, key: "model_level_db", minimum: -60, maximum: 0,
    defaultValue: -20,
  },
  {
    index: 1, key: "field_turbulence", minimum: 0, maximum: 1,
    defaultValue: .65,
  },
];
const event = {
  strength: .8, location: .8, hardness: .65, implement: .75,
  contactSpread: .2, constraint: 0, seed: 17,
};
const analysis = {
  size: 2048, hop: 512, window: "hann", floorDb: -160, dynamicRangeDb: 90,
};
const state = {
  macros: [-36, .8], event, analysis, activeSnapshotId: null,
  reference: { id: "sha256:test", sha256: "test", referenceGainDb: 42 },
};

const fit = snapshotState(state, "Test", descriptors);
assert.equal(fit.schema, "triggerfish.percussion.fit/v1");
assert.equal(fit.renderer.recipe, "metal.cymbal.v1");
assert.equal(fit.renderer.api, 1);
assert.equal(fit.controls.macros, undefined);
assert.deepEqual(fitMacroValues(validateFit(fit, descriptors), descriptors), [-36, .8]);
const roundTrip = validateFit(JSON.parse(JSON.stringify(fit)), descriptors);
const linkedDescriptors = [...descriptors,
  {index:2,key:'bloom_energy_acceleration',minimum:0,maximum:1,defaultValue:.5}];
const splitDescriptors = [...linkedDescriptors,
  {index:3,key:'bloom_energy_sensitivity',minimum:0,maximum:2,defaultValue:1}];
const linked = snapshotState({...state,macros:[-36,.8,.15]},'Linked',linkedDescriptors);
const split = validateFit(linked,splitDescriptors);
assert.deepEqual(fitMacroValues(split,splitDescriptors),[-36,.8,.15,.3]);
assert.equal(linked.instrument.nodes.find(n=>n.id==='body').parameters.bloom_energy_sensitivity,undefined);
split.instrument.nodes.find(n=>n.id==='body').parameters.bloom_energy_sensitivity=.8;
assert.equal(fitMacroValues(validateFit(split,splitDescriptors),splitDescriptors)[3],.8,
  'Existing independent energy sensitivity must not be relinked');
const withBlur = [...descriptors, {index:2,key:'field_phase_tilt',minimum:-2,maximum:2,defaultValue:0}];
const importedBlur = validateFit(fit, withBlur);
assert.deepEqual(fitMacroValues(importedBlur, withBlur), [-36,.8,0]);
assert.equal(fit.instrument.nodes.find(n=>n.id==='body').parameters.field_phase_tilt, undefined,
  'Import does not mutate the saved source');
const oldDoubletDescriptors=[...descriptors,
  {index:2,key:'field_distribution',minimum:0,maximum:3,defaultValue:3},
  {index:3,key:'field_beat_depth',minimum:0,maximum:1,defaultValue:.3},
  {index:4,key:'field_beat_rate_tilt',minimum:-1,maximum:1,defaultValue:.25}];
const oldDoublet=snapshotState({...state,macros:[-36,.8,2,.3,.25]},'Old doublets',oldDoubletDescriptors);
const newDoubletDescriptors=[...oldDoubletDescriptors,{...withBlur[2],index:5}];
assert.deepEqual(fitMacroValues(validateFit(oldDoublet,newDoubletDescriptors),newDoubletDescriptors),
  [-36,.8,2,1,0,0],'Old doublets explicitly retain equal partners and shared gap');
assert.equal(roundTrip.reference.referenceGainDb, 42);
assert.throws(() => validateFit({
  ...fit, reference: { ...fit.reference, referenceGainDb: Infinity },
}, descriptors), /invalid reference gain/);
assert.deepEqual(fitMacroValues(roundTrip, descriptors), [-36, .8]);

const referenceFree = snapshotState(
  { ...state, reference: null }, "No reference", descriptors,
);
assert.equal(validateFit(referenceFree, descriptors).reference, null);

const routedPatch = structuredClone(fit.instrument);
routedPatch.connections[0].enabled = false;
const routedFit = snapshotState({ ...state, patch: routedPatch }, "Routed", descriptors);
assert.equal(validateFit(routedFit, descriptors)
  .instrument.connections[0].enabled, false);

const obsoleteGain = structuredClone(fit);
obsoleteGain.instrument.connections[0].gain = .5;
assert.throws(
  () => validateFit(obsoleteGain, descriptors), /invalid percussion connection/,
);
assert.throws(
  () => validateFit({ ...structuredClone(fit), schema: "old" }, descriptors),
  /unsupported/,
);
const invalidParameter = structuredClone(fit);
invalidParameter.instrument.nodes.find(node => node.id === "body")
  .parameters.field_turbulence = Number.NaN;
assert.throws(
  () => validateFit(invalidParameter, descriptors), /invalid module parameter/,
);
const incompleteEvent = structuredClone(fit);
delete incompleteEvent.controls.event.contactSpread;
assert.throws(() => validateFit(incompleteEvent, descriptors), /invalid event/);

const file = {
  async text() { return JSON.stringify(fit); },
};
assert.deepEqual(
  fitMacroValues(await readFit(file, descriptors), descriptors), [-36, .8],
);

const snareDescriptors = [
  { index: 0, key: "model_level_db", minimum: -60, maximum: 0,
    defaultValue: -14 },
  { index: 1, key: "contact_direct_level", minimum: 0, maximum: 4,
    defaultValue: .273 },
  { index: 2, key: "contact_body_level", minimum: 0, maximum: 4,
    defaultValue: .78 },
  { index: 3, key: "fm_direct_level", minimum: 0, maximum: 3,
    defaultValue: .004 },
  { index: 4, key: "fm_body_level", minimum: 0, maximum: 3,
    defaultValue: .0256 },
  { index: 5, key: "body_level", minimum: 0, maximum: 3,
    defaultValue: .95 },
  { index: 6, key: "wire_level", minimum: 0, maximum: 6,
    defaultValue: .46 },
];
const snareValues = snareDescriptors.map(item => item.defaultValue);
const snareFit = snapshotState({
  ...state, patch: createSnarePatch(snareDescriptors, snareValues),
  macros: snareValues,
}, "Snare", snareDescriptors);
assert.deepEqual(
  fitMacroValues(validateFit(snareFit, snareDescriptors), snareDescriptors),
  snareValues,
);

const kickDescriptors = [
  { index: 0, key: "model_level_db", minimum: -60, maximum: 0,
    defaultValue: -12 },
  { index: 1, key: "thump_pitch_hz", minimum: 25, maximum: 120,
    defaultValue: 52 },
];
const kickValues = [-9, 58];
const kickFit = snapshotState({
  ...state, macros: kickValues,
  patch: createKickPatch(kickDescriptors, kickValues),
}, "Kick", kickDescriptors);
assert.equal(kickFit.renderer.recipe, "drum.kick.v1");
assert.deepEqual(
  fitMacroValues(validateFit(kickFit, kickDescriptors), kickDescriptors),
  kickValues,
);

console.log("workbench state tests passed");

const motionDescriptors = [...descriptors,
  {index:2,key:'field_motion_depth',minimum:0,maximum:1.5,defaultValue:0},
  {index:3,key:'field_motion_rate',minimum:.1,maximum:200,defaultValue:40},
  {index:4,key:'field_motion_sharing',minimum:0,maximum:1,defaultValue:.5}];
const motionFit=validateFit(fit,motionDescriptors);
assert.deepEqual(fitMacroValues(motionFit,motionDescriptors),[-36,.8,0,40,.5]);
assert.equal(fit.instrument.nodes.find(n=>n.id==='body').parameters.field_motion_depth,undefined);
const partial=structuredClone(motionFit);
delete partial.instrument.nodes.find(n=>n.id==='body').parameters.field_motion_rate;
assert.throws(()=>validateFit(partial,motionDescriptors),/parameter/);
const changedMotion=structuredClone(motionFit);
changedMotion.instrument.nodes.find(n=>n.id==='body').parameters.field_motion_depth=.8;
assert.equal(fitMacroValues(validateFit(changedMotion,motionDescriptors),motionDescriptors)[2],.8);
