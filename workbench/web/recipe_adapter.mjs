import {
  createKickPatch, kickRoutingValues, kickValuesFromPatch,
  patchWithKickValues, validateKickPatch,
} from "./compact_kick_patch.mjs";
import {
  createCrashPatch, macroValuesFromPatch, patchWithMacroValues,
  routingValuesFromPatch, validateCrashAdapterPatch,
} from "./metallic_plate_patch.mjs";

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
]);

export function recipeAdapter(recipe) {
  const result = adapters.get(recipe);
  if (!result) throw new Error(`no compiled adapter for ${recipe}`);
  return result;
}
