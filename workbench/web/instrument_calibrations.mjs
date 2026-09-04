import {
  createAcousticKickPatch, membranePresetValues, patchWithMembraneValues,
} from "./membrane_patch.mjs";
import { metallicCalibrationValues } from "./metallic_calibrations.mjs";

export function calibrationParameterValues(calibration, descriptors) {
  if (calibration.parameter_preset === "acoustic-kick")
    return membranePresetValues("acousticKick", descriptors);
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
  if (calibration.parameter_preset !== "acoustic-kick") return fallbackPatch;
  return patchWithMembraneValues(
    createAcousticKickPatch(descriptors), descriptors, values,
  );
}
