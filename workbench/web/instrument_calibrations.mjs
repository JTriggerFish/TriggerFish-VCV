import { patchWithKickValues } from "./kick_patch.mjs";
import { metallicCalibrationValues } from "./metallic_calibrations.mjs";
import KickCalibration from "./kick_calibration.fit.json" with { type: "json" };

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
  return recipe === "drum.kick.v1" ? kickCalibrationValues(descriptors)
    : descriptors.map(item => item.defaultValue);
}

export function recipeStartingEvent(recipe) {
  return recipe === "drum.kick.v1" ? { ...KickCalibration.controls.event } : null;
}

export function calibrationParameterValues(calibration, descriptors) {
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
  if (calibration.parameter_preset !== "kick") return fallbackPatch;
  return patchWithKickValues(KickCalibration.instrument, descriptors, values);
}
