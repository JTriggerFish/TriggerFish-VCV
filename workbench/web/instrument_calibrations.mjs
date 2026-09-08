import { patchWithKickValues } from "./kick_patch.mjs";
import { metallicCalibrationValues } from "./metallic_calibrations.mjs";
import KickCalibration from "./kick_calibration.fit.json" with { type: "json" };
import { checkedCalibrationValues, referenceCalibration, recipeReferenceCalibration } from "./reference_calibration_library.mjs";
import { recipeAdapter } from "./recipe_adapter.mjs";

export function kickCalibrationValues(descriptors) {
  const values = Object.assign({}, ...KickCalibration.instrument.nodes.map(node => node.parameters));
  if (Object.keys(values).length !== descriptors.length)
    throw new Error("Kick calibration does not match the current parameter surface");
  return descriptors.map(item => {
    const value = values[item.key];
    if (!Number.isFinite(value) || value < item.minimum || value > item.maximum)
      throw new Error(`Invalid kick calibration parameter: ${item.key}`);
    return value;
  });
}

export function recipeStartingValues(recipe, descriptors) {
  const fitted = recipeReferenceCalibration(recipe);
  if (fitted) return checkedCalibrationValues(fitted, descriptors);
  return recipe === "drum.kick.v1" ? kickCalibrationValues(descriptors)
    : descriptors.map(item => item.defaultValue);
}

export function recipeStartingEvent(recipe) {
  const fitted = recipeReferenceCalibration(recipe);
  if (fitted) return { ...fitted.controls.event };
  return recipe === "drum.kick.v1" ? { ...KickCalibration.controls.event } : null;
}

// A fitted gesture may differ from the source sample's metadata. Restore it
// after selecting that reference, whose normal browser action sets the event.
export function calibrationEvent(calibration) {
  const fit = referenceCalibration(calibration.id) ??
    (calibration.parameter_preset === "kick" ? KickCalibration : null);
  return fit ? { ...fit.controls.event } : null;
}

export function calibrationParameterValues(calibration, descriptors) {
  const fitted = referenceCalibration(calibration.id);
  if (fitted) {
    if (fitted.instrument.recipe !== calibration.recipe)
      throw new Error("Calibration recipe differs from the selected target");
    return checkedCalibrationValues(fitted, descriptors);
  }
  if (calibration.parameter_preset === "kick")
    return kickCalibrationValues(descriptors);
  const metallic = metallicCalibrationValues(
    calibration.parameter_preset, descriptors,
  );
  if (metallic) return metallic;
  if (calibration.parameter_preset)
    throw new Error(`Unknown calibration preset: ${calibration.parameter_preset}`);
  return descriptors.map(item => item.defaultValue);
}

export function calibrationPatch(
  calibration, descriptors, values, fallbackPatch,
) {
  const fitted = referenceCalibration(calibration.id);
  if (fitted) {
    if (fitted.instrument.recipe !== calibration.recipe)
      throw new Error("Calibration recipe differs from the selected target");
    return recipeAdapter(fitted.instrument.recipe).withValues(fitted.instrument, descriptors, values);
  }
  if (calibration.parameter_preset !== "kick") return fallbackPatch;
  return patchWithKickValues(KickCalibration.instrument, descriptors, values);
}
