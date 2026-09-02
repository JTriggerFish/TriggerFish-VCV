import assert from "node:assert/strict";

import { validateFit } from "../web/state.mjs";

const descriptors = [
  { key: "level", minimum: -60, maximum: 12 },
  { key: "density", minimum: 0, maximum: 1 },
];
const fit = {
  schema: "triggerfish-percussion-fit-v12",
  renderer: {
    graph: "crash-experimental-v12", api: 12, macros: "crash-macros-v12",
  },
  reference: { id: "sha256:test", sha256: "test" },
  controls: {
    macros: [-36, .8],
    event: {
      strength: .8, location: .8, hardness: .65, implement: .75,
      contactSpread: .2, seed: 17,
    },
    analysis: {
      size: 2048, hop: 512, window: "hann", floorDb: -160, dynamicRangeDb: 90,
    },
  },
};

assert.equal(validateFit(structuredClone(fit), descriptors).renderer.api, 12);
assert.throws(
  () => validateFit({ ...structuredClone(fit), schema: "old" }, descriptors),
  /unsupported/,
);
const wrongCount = structuredClone(fit);
wrongCount.controls.macros.pop();
assert.throws(() => validateFit(wrongCount, descriptors), /unsupported/);
const invalidMacro = structuredClone(fit);
invalidMacro.controls.macros[1] = Number.NaN;
assert.throws(() => validateFit(invalidMacro, descriptors), /invalid saved control/);
const incompleteEvent = structuredClone(fit);
delete incompleteEvent.controls.event.contactSpread;
assert.throws(() => validateFit(incompleteEvent, descriptors), /invalid event/);

const fullDescriptors = Array.from({ length: 167 }, (_, index) => ({
  key: `macro_${index}`, minimum: -1000, maximum: 24000,
  defaultValue: index + 0.25,
}));
for (let index = 119; index < 143; ++index) {
  fullDescriptors[index] = {
    key: `resolved_level_${index - 119}`,
    minimum: -72, maximum: 6, defaultValue: -18 + .25 * (index - 119),
  };
}
for (let index = 87; index < 95; ++index) {
  fullDescriptors[index] = {
    key: `body_decay_active_${index - 87}`,
    minimum: 0, maximum: 1, defaultValue: index === 87 || index === 94 ? 1 : 0,
  };
}
const oldFit = structuredClone(fit);
oldFit.schema = "triggerfish-percussion-fit-v8";
oldFit.renderer = {
  graph: "crash-experimental-v8", api: 8, macros: "crash-macros-v8",
};
oldFit.controls.macros = Array.from({ length: 103 }, (_, index) => index);
const migrated = validateFit(oldFit, fullDescriptors);
assert.equal(migrated.renderer.api, 12);
assert.equal(migrated.controls.macros.length, 167);
assert.equal(migrated.controls.macros[3], fullDescriptors[3].defaultValue);
assert.equal(migrated.controls.macros[14], fullDescriptors[14].defaultValue);
assert.equal(migrated.controls.macros[15], 13);
assert.ok(migrated.controls.macros[119] >= -72 &&
  migrated.controls.macros[119] <= 6);
assert.equal(migrated.controls.macros[142], fullDescriptors[142].defaultValue);
assert.equal(migrated.controls.macros[166], fullDescriptors[166].defaultValue);

const v7 = structuredClone(fit);
v7.schema = "triggerfish-percussion-fit-v7";
v7.renderer = {
  graph: "crash-experimental-v7", api: 7, macros: "crash-macros-v7",
};
v7.controls.macros = Array.from({ length: 97 }, (_, index) => index);
const migratedV7 = validateFit(v7, fullDescriptors);
assert.equal(migratedV7.renderer.api, 12);
assert.equal(migratedV7.controls.macros[4], 3);
assert.equal(migratedV7.controls.macros[3], fullDescriptors[3].defaultValue);
assert.equal(migratedV7.controls.macros[7], 5);
assert.equal(migratedV7.controls.macros[8], 6);
assert.equal(migratedV7.controls.macros[15], 7);
assert.equal(migratedV7.controls.macros[95], 73);

const v10 = structuredClone(fit);
v10.schema = "triggerfish-percussion-fit-v10";
v10.renderer = {
  graph: "crash-experimental-v10", api: 10, macros: "crash-macros-v10",
};
// Versions through v11 contained 152 values.  Use a genuinely old-sized
// payload here; feeding a current-sized array would correctly be rejected.
v10.controls.macros = Array(152).fill(.5);
const relative = Array.from({ length: 24 }, (_, index) => index - 11.5);
relative.forEach((value, index) => { v10.controls.macros[104 + index] = value; });
const migratedV10 = validateFit(v10, fullDescriptors);
assert.equal(migratedV10.renderer.api, 12);
assert.equal(migratedV10.controls.macros[119],
  fullDescriptors[119].defaultValue + relative[0]);
assert.equal(migratedV10.controls.macros[142],
  fullDescriptors[142].defaultValue + relative[23]);

const v11 = structuredClone(fit);
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
assert.equal(migratedV11.controls.macros[3], fullDescriptors[3].defaultValue);
assert.equal(migratedV11.controls.macros[4], .3);
assert.equal(migratedV11.controls.macros[72], 150);
assert.equal(migratedV11.controls.macros[79], 12);
assert.equal(migratedV11.controls.macros[80], 12);
assert.deepEqual(migratedV11.controls.macros.slice(87, 95),
  [1, 1, 1, 1, 1, 1, 0, 1]);
assert.equal(migratedV11.controls.macros[95], 440);
assert.equal(migratedV11.controls.macros[143], .7);

console.log("workbench state tests passed");
