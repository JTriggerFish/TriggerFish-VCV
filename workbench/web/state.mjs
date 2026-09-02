const FitSchema = "triggerfish-percussion-fit-v1";

export function snapshotState(state, name = "Snapshot") {
  return Object.freeze({
    schema: FitSchema,
    id: crypto.randomUUID(),
    parentId: state.activeSnapshotId ?? null,
    createdAt: new Date().toISOString(),
    name,
    renderer: { graph: "crash-experimental-v1", api: 1, macros: "crash-macros-v1" },
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

export function validateFit(value) {
  if (value?.schema !== FitSchema || value?.renderer?.api !== 1 ||
      !Array.isArray(value?.controls?.macros)) {
    throw new Error("unsupported or incomplete percussion fit");
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

export async function readFit(file) {
  return validateFit(JSON.parse(await file.text()));
}
