import assert from "node:assert/strict";

import { validateFit } from "../web/state.mjs";

const descriptors = [
  { key: "level", minimum: -60, maximum: 12 },
  { key: "density", minimum: 0, maximum: 1 },
];
const fit = {
  schema: "triggerfish-percussion-fit-v4",
  renderer: {
    graph: "crash-experimental-v4", api: 4, macros: "crash-macros-v4",
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

assert.equal(validateFit(structuredClone(fit), descriptors).renderer.api, 4);
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

console.log("workbench state tests passed");
