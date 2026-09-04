import {
  createAcousticKickPatch, membranePresetValues, patchWithMembraneValues,
} from "./membrane_patch.mjs";

export function calibrationParameterValues(calibration, descriptors) {
  if (calibration.parameter_preset === "acoustic-kick")
    return membranePresetValues("acousticKick", descriptors);
  if (calibration.parameter_preset)
    throw new Error(`Unknown calibration preset: ${calibration.parameter_preset}`);
  return descriptors.map(item => item.defaultValue);
}

export function calibrationPatch(
  calibration, descriptors, values, fallbackPatch,
) {
  if (calibration.parameter_preset !== "acoustic-kick") return fallbackPatch;
  return patchWithMembraneValues(
    createAcousticKickPatch(descriptors), descriptors, values,
  );
}
