const FitSchema = "triggerfish-percussion-fit-v12";
const CommittedFitSchema = "triggerfish-percussion-fit-v7";
const PreviousFitSchema = "triggerfish-percussion-fit-v8";
const IntermediateFitSchema = "triggerfish-percussion-fit-v9";
const PriorFitSchema = "triggerfish-percussion-fit-v10";
const UnifiedFitSchema = "triggerfish-percussion-fit-v11";

export function snapshotState(state, name = "Snapshot") {
  return Object.freeze({
    schema: FitSchema,
    id: crypto.randomUUID(),
    parentId: state.activeSnapshotId ?? null,
    createdAt: new Date().toISOString(),
    name,
    renderer: { graph: "crash-experimental-v12", api: 12, macros: "crash-macros-v12" },
    reference: state.reference ? {
      id: state.reference.id,
      sha256: state.reference.sha256,
      name: state.reference.name,
      sampleRate: state.reference.sampleRate,
      channels: state.reference.channels,
      bits: state.reference.bits,
      duration: state.reference.duration,
      corpus: state.reference.corpus ?? null,
      cell: state.reference.cell ?? null,
    } : null,
    controls: {
      macros: [...state.macros],
      event: { ...state.event },
      analysis: { ...state.analysis },
    },
  });
}

export function validateFit(value, descriptors = []) {
  value = migrateLegacyFit(value, descriptors);
  const macros = value?.controls?.macros;
  const event = value?.controls?.event;
  const analysis = value?.controls?.analysis;
  if (value?.schema !== FitSchema || value?.renderer?.api !== 12 ||
      value?.renderer?.graph !== "crash-experimental-v12" ||
      value?.renderer?.macros !== "crash-macros-v12" ||
      !Array.isArray(macros) || macros.length !== descriptors.length ||
      !event || !analysis || typeof value?.reference?.id !== "string" ||
      typeof value?.reference?.sha256 !== "string") {
    throw new Error("unsupported or incomplete percussion fit");
  }
  descriptors.forEach((descriptor, index) => {
    const macro = macros[index];
    if (!Number.isFinite(macro) || macro < descriptor.minimum ||
        macro > descriptor.maximum) {
      throw new Error(`invalid saved control: ${descriptor.key}`);
    }
  });
  for (const key of [
    "strength", "location", "hardness", "implement", "contactSpread",
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

function migrateLegacyFit(value, descriptors) {
  if (value?.schema === FitSchema) return value;
  if (descriptors.length !== 167) return value;
  const oldDefaults = v11Defaults(descriptors);
  let old;
  if (value?.schema === UnifiedFitSchema && value?.renderer?.api === 11 &&
      Array.isArray(value?.controls?.macros) &&
      value.controls.macros.length === 152) {
    old = [...value.controls.macros];
  }
  if (value?.schema === PriorFitSchema && value?.renderer?.api === 10 &&
      Array.isArray(value?.controls?.macros) && value.controls.macros.length === 152) {
    old = [...value.controls.macros];
    migrateRelativeModeLevels(old, oldDefaults, 24);
  }
  if (value?.schema === IntermediateFitSchema && value?.renderer?.api === 9 &&
      Array.isArray(value?.controls?.macros) && value.controls.macros.length === 128) {
    old = [
      ...value.controls.macros,
      ...oldDefaults.slice(128),
    ];
    migrateRelativeModeLevels(old, oldDefaults, 24);
  }
  if (value?.schema === CommittedFitSchema && value?.renderer?.api === 7 &&
      Array.isArray(value?.controls?.macros) && value.controls.macros.length === 97) {
    const source = value.controls.macros;
    old = [...oldDefaults];
    for (let index = 0; index < 5; ++index) old[index] = source[index];
    old[6] = source[5];
    old[7] = source[6];
    for (let index = 7; index < 73; ++index) old[index + 7] = source[index];
    for (let index = 0; index < 12; ++index) {
      old[80 + index] = source[73 + index];
      old[104 + index] = source[85 + index];
    }
    migrateRelativeModeLevels(old, oldDefaults, 12);
  }
  if (value?.schema === PreviousFitSchema && value?.renderer?.api === 8 &&
      Array.isArray(value?.controls?.macros) && value.controls.macros.length === 103) {
    const source = value.controls.macros;
    old = [...oldDefaults];
    for (let index = 0; index < 13; ++index) old[index] = source[index];
    for (let index = 13; index < 79; ++index) old[index + 1] = source[index];
    for (let index = 0; index < 12; ++index) {
      old[80 + index] = source[79 + index];
      old[104 + index] = source[91 + index];
    }
    migrateRelativeModeLevels(old, oldDefaults, 12);
  }
  return old ? migrated(value, expandV11(old, descriptors)) : value;
}

function v11Defaults(descriptors) {
  const current = descriptors.map(descriptor => descriptor.defaultValue);
  const old = Array(152);
  for (let index = 0; index < 3; ++index) old[index] = current[index];
  old[3] = current[4];
  for (let index = 4; index < 70; ++index) old[index] = current[index + 1];
  const slots = [0, 1, 2, 3, 7];
  slots.forEach((slot, index) => {
    old[70 + index] = current[71 + slot];
    old[75 + index] = current[79 + slot];
  });
  for (let index = 80; index < 152; ++index) old[index] = current[index + 15];
  return old;
}

function expandV11(old, descriptors) {
  const macros = descriptors.map(descriptor => descriptor.defaultValue);
  for (let index = 0; index < 3; ++index) macros[index] = old[index];
  macros[4] = old[3];
  for (let index = 4; index < 70; ++index) macros[index + 1] = old[index];
  // Preserve every old knot as an interior point. The new boundaries repeat
  // the first/last T60, exactly matching the old constant extrapolation.
  macros[79] = old[75];
  macros[87] = 1;
  for (let index = 0; index < 5; ++index) {
    const slot = index + 1;
    macros[71 + slot] = old[70 + index];
    macros[79 + slot] = old[75 + index];
    macros[87 + slot] = 1;
  }
  macros[86] = old[79];
  macros[94] = 1;
  for (let index = 80; index < 152; ++index) macros[index + 15] = old[index];
  return macros;
}

function migrateRelativeModeLevels(macros, defaults, count) {
  const first = 104;
  const mean = macros.slice(first, first + count)
    .reduce((sum, value) => sum + value, 0) / count;
  for (let index = 0; index < count; ++index) {
    macros[first + index] = Math.max(-72, Math.min(
      6, defaults[first + index] + macros[first + index] - mean,
    ));
  }
}

function migrated(value, macros) {
  return {
    ...value,
    schema: FitSchema,
    renderer: {
      graph: "crash-experimental-v12", api: 12, macros: "crash-macros-v12",
    },
    controls: { ...value.controls, macros },
  };
}

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
  return validateFit(JSON.parse(await file.text()), descriptors);
}
