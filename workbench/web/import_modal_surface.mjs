// Explicit, sound-preserving imports of the preceding noisiness surface.
// Nonzero legacy percentage drift cannot be represented by absolute-Hz drift.
export function importModalSurface(value, descriptors) {
  if (value?.instrument?.recipe !== "metal.cymbal.v1" ||
      !descriptors.some(d => d.key === "field_wander_hz")) return value;
  const source = value.instrument.nodes?.find(n => n.id === "body")?.parameters;
  if (!source || (!Object.hasOwn(source, "field_turbulence_centre") &&
                  !Object.hasOwn(source, "field_drift_depth"))) return value;
  const copy = structuredClone(value);
  const p = copy.instrument.nodes.find(n => n.id === "body").parameters;
  if (Object.hasOwn(p, "field_turbulence_centre")) {
    const c=p.field_turbulence_centre, g=p.field_turbulence, s=p.field_turbulence_slope;
    if (![c,g,s].every(Number.isFinite) || c<1 || c>15000 || g<0 || g>4 || Math.abs(s)>1)
      throw Error("Invalid legacy noisiness curve; nothing was clamped or discarded.");
    p.field_turbulence = g * (1000/c)**s;
    delete p.field_turbulence_centre;
  }
  if (Object.hasOwn(p, "field_drift_depth")) {
    if (p.field_drift_depth !== 0)
      throw Error("This fit uses old percentage drift. Its sound cannot be converted exactly to Hz wander. Keep the original file; set old drift depth to zero before importing, then adjust Pitch wander here.");
    if (Object.hasOwn(p, "field_wander_hz") || Object.hasOwn(p, "field_wander_rate"))
      throw Error("Fit mixes old drift and new wander controls.");
    if (!Number.isFinite(p.field_drift_rate) || p.field_drift_rate<.1 || p.field_drift_rate>40)
      throw Error("Invalid legacy drift speed.");
    p.field_wander_hz = 0;
    p.field_wander_rate = p.field_drift_rate;
    delete p.field_drift_depth;
    delete p.field_drift_rate;
  }
  return copy;
}

// The old body EQ becomes the shared final EQ. This is an intentional topology
// conversion, not an exact emulation of two independently filtered sources.
export function importOutputEq(value, descriptors) {
  if (value?.instrument?.recipe !== "metal.cymbal.v1" ||
      !descriptors.some(d => d.key === "output_eq_enabled")) return value;
  const source = value.instrument.nodes?.find(n => n.id === "observation")?.parameters;
  const suffixes = ["radiation_enabled", "low_cut", "colour_frequency", "colour_gain", "high_cut"];
  const oldKeys = ["direct", "body"].flatMap(prefix => suffixes.map(s => `${prefix}_${s}`));
  if (!source || !oldKeys.some(key => Object.hasOwn(source, key))) return value;
  const newKeys = suffixes.map(s => s === "radiation_enabled" ? "output_eq_enabled" : `output_${s}`);
  if (newKeys.some(key => Object.hasOwn(source, key)))
    throw Error("Fit mixes per-path EQ controls with the shared output EQ.");
  const bounds = [[0,1], [10,1000], [100,18000], [-18,18], [1000,22000]];
  for (const [i, key] of oldKeys.entries()) {
    const v = source[key], [low, high] = bounds[i % suffixes.length];
    if (!Number.isFinite(v) || v < low || v > high)
      throw Error(`Missing or invalid old EQ control: ${key}`);
  }
  const copy = structuredClone(value);
  const p = copy.instrument.nodes.find(n => n.id === "observation").parameters;
  suffixes.forEach((s, i) => { p[newKeys[i]] = p[`body_${s}`]; });
  oldKeys.forEach(key => { delete p[key]; });
  return copy;
}
