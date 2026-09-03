import {
  createKickPatch, kickRoutingValues, kickValuesFromPatch,
  patchWithKickValues, validateKickPatch,
} from "./compact_kick_patch.mjs";
import {
  createCrashPatch, macroValuesFromPatch, patchWithMacroValues,
  routingValuesFromPatch, validateCrashAdapterPatch,
} from "./metallic_plate_patch.mjs";
import {
  createMembranePatch, membraneRoutingValues, membraneValuesFromPatch,
  patchWithMembraneValues, validateMembranePatch,
} from "./membrane_patch.mjs";
import {
  createSnarePatch, patchWithSnareValues, snareRoutingValues,
  snareValuesFromPatch, validateSnarePatch,
} from "./snare_patch.mjs";

const adapters = new Map([
  ["metal.cymbal.v1", {
    create: createCrashPatch,
    withValues: patchWithMacroValues,
    values: macroValuesFromPatch,
    routing: routingValuesFromPatch,
    validate: validateCrashAdapterPatch,
  }],
  ["drum.kick-fm.v1", {
    create: createKickPatch,
    withValues: patchWithKickValues,
    values: kickValuesFromPatch,
    routing: kickRoutingValues,
    validate: validateKickPatch,
  }],
  ["drum.membrane.v1", {
    create: createMembranePatch,
    withValues: patchWithMembraneValues,
    values: membraneValuesFromPatch,
    routing: membraneRoutingValues,
    validate: validateMembranePatch,
  }],
  ["drum.snare.v1", {
    create: createSnarePatch,
    withValues: patchWithSnareValues,
    values: snareValuesFromPatch,
    routing: snareRoutingValues,
    validate: validateSnarePatch,
  }],
]);

export function recipeAdapter(recipe) {
  const result = adapters.get(recipe);
  if (!result) throw new Error(`no compiled adapter for ${recipe}`);
  return result;
}
