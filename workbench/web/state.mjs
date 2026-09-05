import { createCrashPatch } from "./metallic_plate_patch.mjs";
import { validatePatch } from "./percussion_patch.mjs";
import { recipeAdapter } from "./recipe_adapter.mjs";

const FitSchema = "triggerfish.percussion.fit/v1";
const RendererApi = 1;
const RendererAdapter = "percussion-recipe-v1";

export function snapshotState(state, name = "Snapshot", descriptors = []) {
  const source = state.patch ?? createCrashPatch(descriptors, state.macros);
  const instrument = recipeAdapter(source.recipe).withValues(
    source, descriptors, state.macros,
  );
  return Object.freeze({
    schema: FitSchema,
    id: crypto.randomUUID(),
    parentId: state.activeSnapshotId ?? null,
    createdAt: new Date().toISOString(),
    name,
    renderer: {
      recipe: instrument.recipe, api: RendererApi, adapter: RendererAdapter,
    },
    instrument,
    reference: state.reference ? {
      id: state.reference.id,
      sha256: state.reference.sha256,
      name: state.reference.name,
      sampleRate: state.reference.sampleRate,
      channels: state.reference.channels,
      bits: state.reference.bits,
      duration: state.reference.duration,
      referenceGainDb: state.reference.referenceGainDb ?? 0,
      corpus: state.reference.corpus ?? null,
      cell: state.reference.cell ?? null,
    } : null,
    controls: {
      event: { ...state.event },
      analysis: { ...state.analysis },
    },
  });
}

export function validateFit(value, descriptors = []) {
  const event = value?.controls?.event;
  const analysis = value?.controls?.analysis;
  const validReference = value?.reference === null ||
    typeof value?.reference?.id === "string" &&
    typeof value?.reference?.sha256 === "string";
  if (value?.schema !== FitSchema || value?.renderer?.api !== RendererApi ||
      value?.renderer?.recipe !== value?.instrument?.recipe ||
      value?.renderer?.adapter !== RendererAdapter ||
      !event || !analysis || !validReference) {
    throw new Error("unsupported or incomplete percussion fit");
  }
  validatePatch(value.instrument, descriptors);
  const referenceGain = value.reference?.referenceGainDb ?? 0;
  if (!Number.isFinite(referenceGain) || Math.abs(referenceGain) > 120)
    throw new Error("invalid reference gain");
  recipeAdapter(value.instrument.recipe).validate(value.instrument);
  for (const key of [
    "strength", "location", "hardness", "implement", "contactSpread",
    "constraint",
  ]) {
    if (!Number.isFinite(event[key]) || event[key] < 0 || event[key] > 1) {
      throw new Error(`invalid event ${key}`);
    }
  }
  if (!Number.isInteger(event.seed) || event.seed < 0 || event.seed > 0xffffffff) {
    throw new Error("invalid event seed");
  }
  if (!Number.isInteger(analysis.size) || analysis.size < 2 ||
      (analysis.size & (analysis.size - 1)) !== 0 ||
      !Number.isInteger(analysis.hop) || analysis.hop < 1 ||
      analysis.hop > analysis.size || !Number.isFinite(analysis.floorDb) ||
      !Number.isFinite(analysis.dynamicRangeDb) ||
      !["hann", "blackman-harris", "rectangular"].includes(analysis.window)) {
    throw new Error("invalid analysis settings");
  }
  return value;
}

export function fitParameterValues(fit, descriptors) {
  return recipeAdapter(fit.instrument.recipe).values(
    fit.instrument, descriptors,
  );
}

export const fitMacroValues = fitParameterValues;

export function downloadFit(snapshot) {
  const blob = new Blob([JSON.stringify(snapshot, null, 2) + "\n"], {
    type: "application/json",
  });
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = `${snapshot.name.replace(/[^a-z0-9_-]+/gi, "-")}.json`;
  link.click();
  setTimeout(() => URL.revokeObjectURL(link.href), 0);
}

export async function readFit(file, descriptors) {
  const value = JSON.parse(await file.text());
  const recipe = value?.instrument?.recipe;
  const resolved = typeof descriptors === "function"
    ? descriptors(recipe) : descriptors;
  return validateFit(value, resolved);
}
