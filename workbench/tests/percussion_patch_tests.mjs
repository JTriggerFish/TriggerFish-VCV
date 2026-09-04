import assert from "node:assert/strict";

import {
  createKickPatch, kickRoutingValues, validateKickPatch,
} from "../web/compact_kick_patch.mjs";
import {
  createCrashPatch, macroValuesFromPatch, routingValuesFromPatch,
  validateCrashAdapterPatch,
} from "../web/metallic_plate_patch.mjs";
import { PatchSchema, validatePatch } from "../web/percussion_patch.mjs";
import {
  createAcousticKickPatch, createMembranePatch, membranePresetValues,
  membraneRoutingValues,
  validateMembranePatch,
} from "../web/membrane_patch.mjs";
import {
  createSnarePatch, snareRoutingValues, validateSnarePatch,
} from "../web/snare_patch.mjs";

const descriptors = [
  {
    index: 0, key: "impact_width", minimum: .25, maximum: 4,
    defaultValue: .65,
  },
  {
    index: 1, key: "bloom_rate", minimum: 0, maximum: 16,
    defaultValue: 2,
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
assert.deepEqual(routingValuesFromPatch(patch), [true, true, true]);
const withoutLayout = structuredClone(patch);
delete withoutLayout.nodes[0].editor;
assert.equal(validatePatch(withoutLayout, descriptors), withoutLayout);

const cyclic = structuredClone(patch);
cyclic.connections.push({
  id: "cycle", from: "output.audio", to: "body.primary", enabled: true,
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
unavailableRecipe.recipe = "metal.pair.v1";
assert.throws(() => validatePatch(unavailableRecipe, descriptors), /unsupported/);

const unsupported = structuredClone(patch);
unsupported.connections.pop();
assert.throws(() => validateCrashAdapterPatch(unsupported), /metallic-plate/);

for (let mask = 0; mask < 8; ++mask) {
  const candidate = structuredClone(patch);
  for (let index = 0; index < 3; ++index) {
    candidate.connections[index].enabled = Boolean(mask & (1 << index));
  }
  const direct = Boolean(mask & (1 << 1));
  const body = Boolean(mask & (1 << 0)) && Boolean(mask & (1 << 2));
  validatePatch(candidate, descriptors);
  if (direct || body) {
    assert.equal(validateCrashAdapterPatch(candidate), candidate);
  } else {
    assert.throws(
      () => validateCrashAdapterPatch(candidate), /no audible route/,
      `route mask ${mask} should be silent`,
    );
  }
  assert.deepEqual(
    routingValuesFromPatch(candidate),
    Array.from({ length: 3 }, (_, index) => Boolean(mask & (1 << index))),
  );
}

const invalidGain = structuredClone(patch);
invalidGain.connections[0].gain = 1;
assert.throws(() => validatePatch(invalidGain, descriptors), /connection/);

const disabledCycle = structuredClone(patch);
disabledCycle.connections.push({
  id: "disabled-cycle", from: "output.audio", to: "body.primary",
  enabled: false,
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
assert.deepEqual(kickRoutingValues(kick), [true, true, true]);
for (let mask = 0; mask < 8; ++mask) {
  const candidate = structuredClone(kick);
  for (let index = 0; index < 3; ++index)
    candidate.connections[index].enabled = Boolean(mask & (1 << index));
  validatePatch(candidate, kickDescriptors);
  if (mask) assert.equal(validateKickPatch(candidate), candidate);
  else assert.throws(() => validateKickPatch(candidate), /no audible route/);
  assert.deepEqual(
    kickRoutingValues(candidate),
    Array.from({ length: 3 }, (_, index) => Boolean(mask & (1 << index))),
  );
}

const brokenKickOutput = structuredClone(kick);
brokenKickOutput.connections[3].enabled = false;
assert.throws(() => validateKickPatch(brokenKickOutput), /required/);

const membraneDescriptors = [
  { index: 0, key: "model_level_db", minimum: -60, maximum: 12, defaultValue: -10 },
  { index: 1, key: "fundamental_hz", minimum: 25, maximum: 500, defaultValue: 105 },
  { index: 2, key: "tension_octaves", minimum: -.25, maximum: .6, defaultValue: .11 },
  { index: 3, key: "contact_direct_level", minimum: 0, maximum: 4, defaultValue: .245 },
  { index: 4, key: "contact_body_level", minimum: 0, maximum: 4, defaultValue: .7 },
  { index: 5, key: "fm_direct_level", minimum: 0, maximum: 3, defaultValue: .0144 },
  { index: 6, key: "fm_body_level", minimum: 0, maximum: 3, defaultValue: .081 },
  { index: 7, key: "direct_level", minimum: 0, maximum: 3, defaultValue: .3 },
  { index: 8, key: "equalizer_mode", minimum: 0, maximum: 2, defaultValue: 1,
    scale: "choice" },
  { index: 9, key: "colour_frequency_hz", minimum: 40, maximum: 20000,
    defaultValue: 2800, scale: "logarithmic" },
  { index: 10, key: "colour_gain_db", minimum: -24, maximum: 24,
    defaultValue: 0 },
];
const membrane = createMembranePatch(
  membraneDescriptors,
  [-10, 105, .11, .245, .7, .0144, .081, .3, 1, 2800, 0],
);
assert.equal(validatePatch(membrane, membraneDescriptors), membrane);
assert.equal(validateMembranePatch(membrane), membrane);
assert.deepEqual(membraneRoutingValues(membrane),
  [true, true, true, true, true]);
const acousticKick = membranePresetValues(
  "acousticKick", membraneDescriptors,
);
assert.equal(acousticKick[1], 35);
assert.equal(acousticKick[3], 1.64);
assert.equal(acousticKick[6], .05625);
const acousticKickPatch = createAcousticKickPatch(
  membraneDescriptors,
);
assert.equal(acousticKickPatch.id, "factory.membrane.acoustic-kick-01");
assert.equal(acousticKickPatch.name, "Acoustic kick");
assert.equal(validateMembranePatch(acousticKickPatch), acousticKickPatch);
const tomAfterKick = membranePresetValues("tom", membraneDescriptors);
assert.deepEqual(acousticKick.slice(9), [600, 10]);
assert.deepEqual(tomAfterKick.slice(9), [2800, 0],
  "factory presets restore their own defaults instead of prior values");
const fractionalChoice = structuredClone(membrane);
fractionalChoice.nodes.find(node => node.id === "membrane-eq")
  .parameters.equalizer_mode = 1.5;
assert.throws(
  () => validatePatch(fractionalChoice, membraneDescriptors), /parameter/,
);
const renamedRoutes = structuredClone(membrane);
renamedRoutes.connections.slice(0, 5).forEach((route, index) => {
  route.id = `arbitrary-${index}`;
});
renamedRoutes.connections[0].enabled = false;
assert.equal(validateMembranePatch(renamedRoutes), renamedRoutes);
assert.deepEqual(membraneRoutingValues(renamedRoutes),
  [false, true, true, true, true]);
const brokenMembrane = structuredClone(membrane);
brokenMembrane.connections.find(route => route.required).enabled = false;
assert.throws(() => validateMembranePatch(brokenMembrane), /required/);
for (let mask = 0; mask < 32; ++mask) {
  const candidate = structuredClone(membrane);
  for (let index = 0; index < 5; ++index)
    candidate.connections[index].enabled = Boolean(mask & (1 << index));
  validatePatch(candidate, membraneDescriptors);
  const direct = Boolean(mask & 1) || Boolean(mask & 4);
  const body = Boolean(mask & 16) &&
    (Boolean(mask & 2) || Boolean(mask & 8));
  if (direct || body) assert.equal(validateMembranePatch(candidate), candidate);
  else assert.throws(
    () => validateMembranePatch(candidate), /no audible route/,
  );
}

const snareDescriptors = [
  ...membraneDescriptors,
  { index: 11, key: "wire_level", minimum: 0, maximum: 3, defaultValue: .46 },
  { index: 12, key: "wire_sensitivity", minimum: 0, maximum: 32, defaultValue: 9 },
  { index: 13, key: "wire_density", minimum: 0, maximum: 1, defaultValue: .9 },
];
const snare = createSnarePatch(
  snareDescriptors,
  [-10, 185, .08, .273, .78, .004, .0256, .24, 1, 3200, 1.5,
    .46, 9, .9],
);
assert.equal(validatePatch(snare, snareDescriptors), snare);
assert.equal(validateSnarePatch(snare), snare);
assert.deepEqual(snareRoutingValues(snare),
  [true, true, true, true, true, true, true]);
for (let mask = 0; mask < 128; ++mask) {
  const candidate = structuredClone(snare);
  for (let index = 0; index < 7; ++index)
    candidate.connections[index].enabled = Boolean(mask & (1 << index));
  validatePatch(candidate, snareDescriptors);
  const direct = Boolean(mask & 1) || Boolean(mask & 4);
  const bodyDrive = Boolean(mask & 2) || Boolean(mask & 8);
  const body = bodyDrive && Boolean(mask & 16);
  const wires = bodyDrive && Boolean(mask & 32) && Boolean(mask & 64);
  if (direct || body || wires)
    assert.equal(validateSnarePatch(candidate), candidate);
  else assert.throws(() => validateSnarePatch(candidate), /no audible route/);
}
const brokenSnare = structuredClone(snare);
brokenSnare.connections.find(route => route.required).enabled = false;
assert.throws(() => validateSnarePatch(brokenSnare), /required/);

console.log("percussion patch tests passed");
