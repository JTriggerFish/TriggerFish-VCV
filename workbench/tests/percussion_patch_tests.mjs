import assert from "node:assert/strict";

import {
  createKickPatch, kickRoutingValues, validateKickPatch,
} from "../web/compact_kick_patch.mjs";
import {
  createCrashPatch, macroValuesFromPatch, routingValuesFromPatch,
  validateCrashAdapterPatch,
} from "../web/metallic_plate_patch.mjs";
import { PatchSchema, validatePatch } from "../web/percussion_patch.mjs";

const descriptors = [
  {
    index: 0, key: "impact_width", minimum: .25, maximum: 4,
    defaultValue: .65,
  },
  {
    index: 1, key: "bloom_level", minimum: 0, maximum: 2,
    defaultValue: 1,
  },
  {
    index: 2, key: "field_turbulence", minimum: 0, maximum: 1,
    defaultValue: .5,
  },
  {
    index: 3, key: "model_level_db", minimum: -60, maximum: 12,
    defaultValue: -20,
  },
];

const patch = createCrashPatch(descriptors, [.8, 1.2, .7, -18]);
assert.equal(patch.schema, PatchSchema);
assert.equal(validatePatch(patch, descriptors), patch);
assert.equal(validateCrashAdapterPatch(patch), patch);
assert.deepEqual(patch.performanceControls,
  ["strength", "location", "hardness", "implement", "contactSpread"]);
assert.deepEqual(macroValuesFromPatch(patch, descriptors), [.8, 1.2, .7, -18]);
assert.ok(patch.nodes.every(node => !Array.isArray(node.parameters)));
assert.deepEqual(routingValuesFromPatch(patch), [1, 1, 1, 1, 1]);
const withoutLayout = structuredClone(patch);
delete withoutLayout.nodes[0].editor;
assert.equal(validatePatch(withoutLayout, descriptors), withoutLayout);

const cyclic = structuredClone(patch);
cyclic.connections.push({
  id: "cycle", from: "output.audio", to: "bloom.drive",
});
assert.throws(() => validatePatch(cyclic, descriptors), /cycle/);

const badPort = structuredClone(patch);
badPort.connections[0].to = "body.missing";
assert.throws(() => validatePatch(badPort, descriptors), /connection/);

const badParameter = structuredClone(patch);
badParameter.nodes.find(node => node.id === "contact")
  .parameters.impact_width = 99;
assert.throws(() => validatePatch(badParameter, descriptors), /parameter/);

const missingParameter = structuredClone(patch);
delete missingParameter.nodes.find(node => node.id === "contact")
  .parameters.impact_width;
assert.throws(() => validatePatch(missingParameter, descriptors), /incomplete/);

const futureEngine = structuredClone(patch);
futureEngine.engineMinimum = 2;
assert.throws(() => validatePatch(futureEngine, descriptors), /unsupported/);

const invalidPerformanceControl = structuredClone(patch);
invalidPerformanceControl.performanceControls.push("unknown");
assert.throws(
  () => validatePatch(invalidPerformanceControl, descriptors), /performance/,
);
const duplicatePerformanceControl = structuredClone(patch);
duplicatePerformanceControl.performanceControls.push("strength");
assert.throws(
  () => validatePatch(duplicatePerformanceControl, descriptors), /performance/,
);

const wrongOwner = structuredClone(patch);
const impactWidth = wrongOwner.nodes.find(node => node.id === "contact")
  .parameters.impact_width;
delete wrongOwner.nodes.find(node => node.id === "contact")
  .parameters.impact_width;
wrongOwner.nodes.find(node => node.id === "body")
  .parameters.impact_width = impactWidth;
assert.equal(validatePatch(wrongOwner, descriptors), wrongOwner);
assert.throws(
  () => validateCrashAdapterPatch(wrongOwner), /parameter owner/,
);

const unavailableRecipe = structuredClone(patch);
unavailableRecipe.recipe = "drum.snare.v1";
assert.throws(() => validatePatch(unavailableRecipe, descriptors), /unsupported/);

const unsupported = structuredClone(patch);
unsupported.connections.pop();
assert.throws(() => validateCrashAdapterPatch(unsupported), /metallic-plate/);

for (let mask = 0; mask < 32; ++mask) {
  const candidate = structuredClone(patch);
  for (let index = 0; index < 5; ++index) {
    candidate.connections[index].enabled = Boolean(mask & (1 << index));
  }
  const direct = Boolean(mask & (1 << 3));
  const immediateBody = Boolean(mask & (1 << 1)) &&
    Boolean(mask & (1 << 4));
  const bloomedBody = Boolean(mask & (1 << 0)) &&
    Boolean(mask & (1 << 2)) && Boolean(mask & (1 << 4));
  validatePatch(candidate, descriptors);
  if (direct || immediateBody || bloomedBody) {
    assert.equal(validateCrashAdapterPatch(candidate), candidate);
  } else {
    assert.throws(
      () => validateCrashAdapterPatch(candidate), /no audible route/,
      `route mask ${mask} should be silent`,
    );
  }
  assert.deepEqual(
    routingValuesFromPatch(candidate),
    Array.from({ length: 5 }, (_, index) => mask & (1 << index) ? 1 : 0),
  );
}

const gained = structuredClone(patch);
gained.connections[2].gain = .375;
assert.equal(validatePatch(gained, descriptors), gained);
assert.equal(validateCrashAdapterPatch(gained), gained);
assert.equal(routingValuesFromPatch(gained)[2], .375);

const invalidGain = structuredClone(patch);
invalidGain.connections[0].gain = 2.1;
assert.throws(() => validatePatch(invalidGain, descriptors), /connection/);

const disabledCycle = structuredClone(patch);
disabledCycle.connections.push({
  id: "disabled-cycle", from: "output.audio", to: "bloom.drive",
  enabled: false, gain: 1,
});
assert.equal(validatePatch(disabledCycle, descriptors), disabledCycle);

const kickDescriptors = [
  {
    index: 0, key: "model_level_db", minimum: -60, maximum: 12,
    defaultValue: -12,
  },
  {
    index: 1, key: "fundamental_hz", minimum: 25, maximum: 120,
    defaultValue: 52,
  },
  {
    index: 2, key: "secondary_ratio", minimum: .5, maximum: 3,
    defaultValue: 1.52,
  },
  {
    index: 3, key: "click_level", minimum: 0, maximum: 1.5,
    defaultValue: .16,
  },
  {
    index: 4, key: "low_cut_hz", minimum: 5, maximum: 500,
    defaultValue: 18,
  },
];
const kick = createKickPatch(kickDescriptors, [-10, 55, 1.6, .2, 20]);
assert.equal(validatePatch(kick, kickDescriptors), kick);
assert.equal(validateKickPatch(kick), kick);
assert.deepEqual(kickRoutingValues(kick), [1, 1, 1]);
for (let mask = 0; mask < 8; ++mask) {
  const candidate = structuredClone(kick);
  for (let index = 0; index < 3; ++index)
    candidate.connections[index].enabled = Boolean(mask & (1 << index));
  validatePatch(candidate, kickDescriptors);
  if (mask) assert.equal(validateKickPatch(candidate), candidate);
  else assert.throws(() => validateKickPatch(candidate), /no audible route/);
  assert.deepEqual(
    kickRoutingValues(candidate),
    Array.from({ length: 3 }, (_, index) => mask & (1 << index) ? 1 : 0),
  );
}

const brokenKickOutput = structuredClone(kick);
brokenKickOutput.connections[3].enabled = false;
assert.throws(() => validateKickPatch(brokenKickOutput), /required/);

console.log("percussion patch tests passed");
