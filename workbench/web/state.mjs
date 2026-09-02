const FitSchema = "triggerfish-percussion-fit-v4";

export function snapshotState(state, name = "Snapshot") {
  return Object.freeze({
    schema: FitSchema,
    id: crypto.randomUUID(),
    parentId: state.activeSnapshotId ?? null,
    createdAt: new Date().toISOString(),
    name,
    renderer: { graph: "crash-experimental-v4", api: 4, macros: "crash-macros-v4" },
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
  const macros = value?.controls?.macros;
  const event = value?.controls?.event;
  const analysis = value?.controls?.analysis;
  if (value?.schema !== FitSchema || value?.renderer?.api !== 4 ||
      value?.renderer?.graph !== "crash-experimental-v4" ||
      value?.renderer?.macros !== "crash-macros-v4" ||
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
