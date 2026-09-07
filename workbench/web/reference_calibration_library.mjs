// Full editable patches, not another layer of overrides or gain correction.
import Snare from "./snare_calibration.fit.json" with {type:"json"};
import HiHat from "./hihat_calibration.fit.json" with {type:"json"};
import Crash from "./crash_calibration.fit.json" with {type:"json"};
import Gong from "./gong_calibration.fit.json" with {type:"json"};
import Ride from "./ride_calibration.fit.json" with {type:"json"};

const Fits = new Map([
  ["snare-standard", Snare], ["hihat-standard", HiHat],
  ["crash-standard", Crash], ["gong-standard", Gong],
  ["ride-standard", Ride],
]);

export function referenceCalibration(id) {
  return Fits.get(id) ?? null;
}

export function recipeReferenceCalibration(recipe) {
  const id = recipe === "drum.snare.v1" ? "snare-standard"
    : recipe === "metal.cymbal.v1" ? "crash-standard" : null;
  return referenceCalibration(id);
}

export function checkedCalibrationValues(fit, descriptors) {
  const values = {};
  for (const node of fit.instrument.nodes) for (const [key, value] of Object.entries(node.parameters)) {
    if (Object.hasOwn(values, key)) throw new Error(`Duplicate calibration parameter: ${key}`);
    values[key] = value;
  }
  if (Object.keys(values).length !== descriptors.length)
    throw new Error("Calibration does not match the current parameter surface");
  return descriptors.map(item => {
    const value = values[item.key];
    if (!Number.isFinite(value) || value < item.minimum || value > item.maximum)
      throw new Error(`Invalid calibration parameter: ${item.key}`);
    return value;
  });
}
