// Design tools only: these write ordinary editable modes, never runtime rules.
// Lowest 16 distinct roots j_(m,n) of J_m, divided by j_(0,1).
// Generated offline with scipy.special.jn_zeros (m=0..11, n=1..8), sorted.
// Angular degeneracy is not duplicated; these are editable design handles.
export const MembraneRatios = Object.freeze([
  1, 1.593340506, 2.135548787, 2.295417267,
  2.653066405, 2.917295455, 3.155464815, 3.500147490,
  3.598484674, 3.647451179, 4.058931883, 4.131738160,
  4.230439128, 4.601044534, 4.610051645, 4.831885263,
]);

export function modalTemplate({
  family = "membrane", fundamental = 55, count = 16,
  level = 0, rolloff = 6, minimumFrequency = 20, maximumFrequency = 15000,
} = {}) {
  if (!["membrane", "harmonic"].includes(family) ||
      ![fundamental, count, level, rolloff, minimumFrequency, maximumFrequency].every(Number.isFinite) ||
      fundamental <= 0 || minimumFrequency <= 0 || maximumFrequency < minimumFrequency ||
      count < 1 || count > 32 || !Number.isInteger(count))
    throw Error("Invalid modal template settings");
  // No invented continuation of the root-ratio list beyond its 16 entries.
  const length = family === "membrane" ? Math.min(count, MembraneRatios.length) : count;
  return Array.from({length}, (_, index) => {
    const ratio = family === "membrane" ? MembraneRatios[index] : index + 1;
    const position = length > 1 ? index / (length - 1) : 0;
    return {
      frequency: fundamental * ratio,
      level: level - rolloff * Math.log2(ratio),
      centre: index === 0 ? 1 : .65 * (index % 2 ? -1 : 1) * (1 - .45 * position),
      edge: .18 + .82 * position,
      turbulence: 0, active: true,
    };
  }).filter(point => point.frequency >= minimumFrequency && point.frequency <= maximumFrequency);
}

export function mountModalTemplates(parent, {capacity, minimumFrequency, maximumFrequency, apply}) {
  parent.replaceChildren();
  const details = document.createElement("details");
  const summary = document.createElement("summary");
  summary.textContent = "Generate editable modes";
  details.append(summary);
  const controls = document.createElement("div");
  controls.className = "modal-toolbar";
  const fields = {};
  for (const [key, label, initial, min, max, step] of [
    ["fundamental", "Fundamental (Hz)", 55, minimumFrequency, maximumFrequency, .1],
    ["count", "Modes", Math.min(16, capacity), 1, capacity, 1],
    ["level", "Top level (dB)", 0, -60, 6, .5],
    ["rolloff", "Falloff (dB/oct)", 6, -12, 24, .5],
  ]) {
    const row = document.createElement("label");
    row.textContent = label + " ";
    const input = document.createElement("input");
    Object.assign(input, {type:"number", value:initial, min, max, step});
    input.style.width = "6em";
    input.ondblclick = () => { input.value = initial; };
    row.append(input); controls.append(row); fields[key] = input;
  }
  const family = document.createElement("select");
  family.setAttribute("aria-label", "Modal formula");
  family.append(new Option("Circular membrane", "membrane"), new Option("Harmonic series", "harmonic"));
  controls.prepend(family);
  const button = document.createElement("button");
  button.textContent = "Replace modes";
  const status = document.createElement("span");
  button.onclick = () => {
    try {
      const values = Object.fromEntries(Object.entries(fields).map(([key, input]) => [key, Number(input.value)]));
      if (Object.values(fields).some(input => !input.checkValidity())) throw Error("Check the input ranges");
      const points = modalTemplate({...values, family:family.value, minimumFrequency, maximumFrequency});
      if (!points.length) throw Error("No modes fall inside the editor frequency range");
      apply(points); status.textContent = points.length + " editable modes written";
    } catch (error) { status.textContent = error.message; }
  };
  controls.append(button, status);
  const note = document.createElement("p");
  note.className = "control-help";
  note.textContent = "Replaces the mode list, not the damping curve. Membrane uses up to 16 root ratios. The formula stops here: every resulting mode stays independently editable and is saved in the patch.";
  details.append(controls, note); parent.append(details);
}
