import assert from "node:assert/strict";

import {
  fitMacroValues, readFit, snapshotState, validateFit,
} from "../web/state.mjs";
import { createKickPatch } from "../web/compact_kick_patch.mjs";

const descriptors = [
  {
    index: 0, key: "model_level_db", minimum: -60, maximum: 12,
    defaultValue: -20,
  },
  {
    index: 1, key: "field_turbulence", minimum: 0, maximum: 1,
    defaultValue: .65,
  },
];
const event = {
  strength: .8, location: .8, hardness: .65, implement: .75,
  contactSpread: .2, seed: 17,
};
const analysis = {
  size: 2048, hop: 512, window: "hann", floorDb: -160, dynamicRangeDb: 90,
};

const state = {
  macros: [-36, .8], event, analysis, activeSnapshotId: null,
  reference: { id: "sha256:test", sha256: "test" },
};
const fit = snapshotState(state, "Test", descriptors);
assert.equal(fit.schema, "triggerfish-percussion-fit-v14");
assert.equal(fit.renderer.recipe, "metal.cymbal.v1");
assert.equal(fit.renderer.api, 14);
assert.equal(fit.controls.macros, undefined);
assert.deepEqual(fitMacroValues(validateFit(fit, descriptors), descriptors), [-36, .8]);
const roundTrip = validateFit(JSON.parse(JSON.stringify(fit)), descriptors);
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

const earlyStructured = structuredClone(fit);
earlyStructured.schema = "triggerfish-percussion-fit-v13";
earlyStructured.renderer.api = 13;
earlyStructured.renderer.adapter = "crash-macros-routes-v13";
earlyStructured.instrument.connections.forEach(connection => {
  delete connection.enabled;
  delete connection.gain;
});
earlyStructured.instrument.nodes.find(node => node.id === "body")
  .parameters.body_tone_wash = .5;
const migratedStructured = validateFit(earlyStructured, descriptors);
assert.equal(migratedStructured.renderer.api, 14);
assert.ok(migratedStructured.instrument.connections.every(connection =>
  connection.enabled && connection.gain === 1));
assert.equal(migratedStructured.instrument.nodes.find(node => node.id === "body")
  .parameters.body_tone_wash, undefined);
const mistypedStructured = structuredClone(earlyStructured);
mistypedStructured.instrument.nodes.find(node => node.id === "body")
  .parameters.body_tone_wahs = .5;
assert.throws(
  () => validateFit(mistypedStructured, descriptors), /module parameter/,
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

const legacyDescriptors = Array.from({ length: 167 }, (_, index) => ({
  index, key: `unused_${index}`, minimum: -1000, maximum: 24000,
  defaultValue: index + .25,
}));
const setDescriptor = (index, key, minimum = -1000, maximum = 24000) => {
  legacyDescriptors[index] = {
    index, key, minimum, maximum,
    defaultValue: legacyDescriptors[index].defaultValue,
  };
};
setDescriptor(0, "model_level_db", -60, 12);
setDescriptor(1, "impact_tone_noise");
setDescriptor(2, "impact_width");
setDescriptor(3, "bloom_level", 0, 24000);
setDescriptor(4, "bloom_nonlinearity", 0, 24000);
setDescriptor(5, "bloom_development");
setDescriptor(6, "bloom_diffusion");
setDescriptor(8, "body_brightness");
setDescriptor(10, "field_turbulence");
setDescriptor(11, "field_packet_spread");
setDescriptor(12, "field_phase_bandwidth");
setDescriptor(13, "field_exchange");
setDescriptor(14, "field_gain");
setDescriptor(15, "direct_gain");
const radiationKeys = [
  "radiation_enabled", "low_cut", "low_cut_q", "colour_frequency",
  "colour_gain", "colour_q", "high_cut", "high_cut_q",
];
radiationKeys.forEach((key, index) => setDescriptor(25 + index, `direct_${key}`));
radiationKeys.forEach((key, index) => setDescriptor(41 + index, `dense_${key}`));
for (let index = 71; index < 79; ++index)
  setDescriptor(index, `body_decay_frequency_${index - 71}`);
for (let index = 79; index < 87; ++index)
  setDescriptor(index, `body_decay_seconds_${index - 79}`);
for (let index = 87; index < 95; ++index) {
  setDescriptor(index, `body_decay_active_${index - 87}`, 0, 1);
  legacyDescriptors[index].defaultValue = index === 87 || index === 94 ? 1 : 0;
}
for (let index = 95; index < 119; ++index)
  setDescriptor(index, `resolved_frequency_${index - 95}`);
for (let index = 119; index < 143; ++index) {
  setDescriptor(index, `resolved_level_${index - 119}`, -72, 6);
  legacyDescriptors[index].defaultValue = -18 + .25 * (index - 119);
}
for (let index = 143; index < 167; ++index)
  setDescriptor(index, `resolved_turbulence_${index - 143}`);

const activeLegacyIndices = [
  ...Array.from({ length: 7 }, (_, index) => index), 8,
  ...Array.from({ length: 6 }, (_, index) => 10 + index),
  ...Array.from({ length: 8 }, (_, index) => 25 + index),
  ...Array.from({ length: 8 }, (_, index) => 41 + index),
  ...Array.from({ length: 96 }, (_, index) => 71 + index),
];
const fullDescriptors = activeLegacyIndices.map((legacyIndex, index) => ({
  ...legacyDescriptors[legacyIndex], index,
}));
const defaults = fullDescriptors.map(item => item.defaultValue);
const valueFor = (values, key) =>
  values[fullDescriptors.findIndex(descriptor => descriptor.key === key)];

const v12Macros = legacyDescriptors.map(item => item.defaultValue);
v12Macros[0] = -31;
v12Macros[10] = .42;
const v12 = {
  schema: "triggerfish-percussion-fit-v12",
  renderer: {
    graph: "crash-experimental-v12", api: 12, macros: "crash-macros-v12",
  },
  reference: { id: "sha256:test", sha256: "test" },
  controls: { macros: v12Macros, event, analysis },
};
const migratedV12 = validateFit(v12, fullDescriptors);
assert.equal(migratedV12.schema, "triggerfish-percussion-fit-v14");
const v12Values = fitMacroValues(migratedV12, fullDescriptors);
assert.equal(valueFor(v12Values, "model_level_db"), -31);
assert.equal(valueFor(v12Values, "field_turbulence"), .42);
let resolvedRecipe;
const loadedV12 = await readFit({ text: async () => JSON.stringify(v12) },
  recipe => {
    resolvedRecipe = recipe;
    return fullDescriptors;
  });
assert.equal(resolvedRecipe, "metal.cymbal.v1");
assert.equal(loadedV12.schema, "triggerfish-percussion-fit-v14");

const oldFit = structuredClone(v12);
oldFit.schema = "triggerfish-percussion-fit-v8";
oldFit.renderer = {
  graph: "crash-experimental-v8", api: 8, macros: "crash-macros-v8",
};
oldFit.controls.macros = Array.from({ length: 103 }, (_, index) => index);
const migrated = validateFit(oldFit, fullDescriptors);
const migratedValues = fitMacroValues(migrated, fullDescriptors);
assert.equal(migrated.renderer.api, 14);
assert.equal(migratedValues.length, 126);
assert.equal(valueFor(migratedValues, "bloom_level"),
  valueFor(defaults, "bloom_level"));
assert.equal(valueFor(migratedValues, "direct_gain"), 13);
assert.ok(valueFor(migratedValues, "resolved_level_0") >= -72 &&
  valueFor(migratedValues, "resolved_level_0") <= 6);

const v7 = structuredClone(v12);
v7.schema = "triggerfish-percussion-fit-v7";
v7.renderer = {
  graph: "crash-experimental-v7", api: 7, macros: "crash-macros-v7",
};
v7.controls.macros = Array.from({ length: 97 }, (_, index) => index);
const v7Values = fitMacroValues(validateFit(v7, fullDescriptors), fullDescriptors);
assert.equal(valueFor(v7Values, "bloom_nonlinearity"), 3);
assert.equal(valueFor(v7Values, "bloom_level"),
  valueFor(defaults, "bloom_level"));
assert.equal(valueFor(v7Values, "body_brightness"), 6);
assert.equal(valueFor(v7Values, "direct_gain"), 7);
assert.equal(valueFor(v7Values, "resolved_frequency_0"), 73);

const v10 = structuredClone(v12);
v10.schema = "triggerfish-percussion-fit-v10";
v10.renderer = {
  graph: "crash-experimental-v10", api: 10, macros: "crash-macros-v10",
};
v10.controls.macros = Array(152).fill(.5);
const relative = Array.from({ length: 24 }, (_, index) => index - 11.5);
relative.forEach((value, index) => { v10.controls.macros[104 + index] = value; });
const v10Values = fitMacroValues(
  validateFit(v10, fullDescriptors), fullDescriptors,
);
assert.equal(valueFor(v10Values, "resolved_level_0"),
  valueFor(defaults, "resolved_level_0") + relative[0]);
assert.equal(valueFor(v10Values, "resolved_level_23"),
  valueFor(defaults, "resolved_level_23") + relative[23]);

const v11 = structuredClone(v12);
v11.schema = "triggerfish-percussion-fit-v11";
v11.renderer = {
  graph: "crash-experimental-v11", api: 11, macros: "crash-macros-v11",
};
v11.controls.macros = Array(152).fill(.5);
v11.controls.macros[3] = .3;
v11.controls.macros[70] = 150;
v11.controls.macros[75] = 12;
v11.controls.macros[80] = 440;
v11.controls.macros[128] = .7;
const migratedV11 = validateFit(v11, fullDescriptors);
const v11Values = fitMacroValues(migratedV11, fullDescriptors);
assert.equal(valueFor(v11Values, "bloom_level"),
  valueFor(defaults, "bloom_level"));
assert.equal(valueFor(v11Values, "bloom_nonlinearity"), .3);
assert.equal(valueFor(v11Values, "body_decay_frequency_1"), 150);
assert.equal(valueFor(v11Values, "body_decay_seconds_0"), 12);
assert.equal(valueFor(v11Values, "body_decay_seconds_1"), 12);
assert.deepEqual(Array.from({ length: 8 }, (_, index) =>
  valueFor(v11Values, `body_decay_active_${index}`)),
  [1, 1, 1, 1, 1, 1, 0, 1]);
assert.equal(valueFor(v11Values, "resolved_frequency_0"), 440);
assert.equal(valueFor(v11Values, "resolved_turbulence_0"), .7);

const kickDescriptors = [
  {
    index: 0, key: "model_level_db", minimum: -60, maximum: 12,
    defaultValue: -12,
  },
  {
    index: 1, key: "fundamental_hz", minimum: 25, maximum: 120,
    defaultValue: 52,
  },
];
const kickState = {
  ...state, macros: [-9, 58],
  patch: createKickPatch(kickDescriptors, [-9, 58]),
};
const kickFit = snapshotState(kickState, "Kick", kickDescriptors);
assert.equal(kickFit.renderer.recipe, "drum.kick-fm.v1");
assert.equal(kickFit.renderer.api, 14);
assert.deepEqual(
  fitMacroValues(validateFit(kickFit, kickDescriptors), kickDescriptors),
  [-9, 58],
);

console.log("workbench state tests passed");
