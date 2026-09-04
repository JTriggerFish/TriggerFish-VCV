export const ModuleRoles = Object.freeze({
  exciter: { name: "Exciter", colour: "#e8b45a" },
  transform: { name: "Transformation", colour: "#b48ae6" },
  body: { name: "Body", colour: "#68a7d8" },
  interaction: { name: "Interaction", colour: "#56d39b" },
  observation: { name: "Observation", colour: "#db7da8" },
  output: { name: "Output", colour: "#9aa6b5" },
});

const types = [
  {
    type: "exciter.contact", version: 1, name: "Contact", role: "exciter",
    inputs: [], outputs: [
      { name: "direct", type: "audio" }, { name: "body", type: "audio" },
      { name: "event", type: "event" },
    ],
  },
  {
    type: "transform.dispersion-loop", version: 1,
    name: "Bloom", role: "transform",
    inputs: [{ name: "drive", type: "audio" }],
    outputs: [{ name: "return", type: "audio" }],
  },
  {
    type: "body.stochastic-modal-field", version: 1,
    name: "Metallic body", role: "body",
    inputs: [{ name: "primary", type: "audio" }],
    outputs: [{ name: "audio", type: "audio" }],
  },
  {
    type: "observation.dual-source", version: 1,
    name: "Observation", role: "observation",
    inputs: [
      { name: "direct", type: "audio" }, { name: "body", type: "audio" },
    ],
    outputs: [{ name: "audio", type: "audio" }],
  },
  {
    type: "output.mono", version: 1, name: "Mono", role: "output",
    inputs: [{ name: "audio", type: "audio" }],
    outputs: [{ name: "audio", type: "audio" }],
  },
  {
    type: "exciter.correlated-fm", version: 1,
    name: "Correlated FM", role: "exciter",
    inputs: [], outputs: [{ name: "audio", type: "audio" }],
  },
  {
    type: "exciter.noise-burst", version: 1,
    name: "Noise burst", role: "exciter",
    inputs: [], outputs: [{ name: "audio", type: "audio" }],
  },
  {
    type: "transform.sum3", version: 1, name: "Source mix", role: "transform",
    inputs: [
      { name: "a", type: "audio" }, { name: "b", type: "audio" },
      { name: "c", type: "audio" },
    ],
    outputs: [{ name: "audio", type: "audio" }],
  },
  {
    type: "observation.single-source", version: 1,
    name: "Observation", role: "observation",
    inputs: [{ name: "source", type: "audio" }],
    outputs: [{ name: "audio", type: "audio" }],
  },
  {
    type: "transform.sum2", version: 1, name: "Two-source mix", role: "transform",
    inputs: [{ name: "a", type: "audio" }, { name: "b", type: "audio" }],
    outputs: [{ name: "audio", type: "audio" }],
  },
  {
    type: "interaction.strike-energy", version: 1,
    name: "Strike energy", role: "interaction",
    inputs: [{ name: "strike", type: "event" }],
    outputs: [{ name: "tension", type: "control" }],
  },
  {
    type: "body.membrane-modal", version: 1,
    name: "Modal membrane", role: "body",
    inputs: [
      { name: "drive", type: "audio" }, { name: "tension", type: "control" },
    ],
    outputs: [{ name: "audio", type: "audio" }],
  },
  {
    type: "observation.equalizer", version: 1,
    name: "Selectable EQ", role: "observation",
    inputs: [{ name: "audio", type: "audio" }],
    outputs: [{ name: "audio", type: "audio" }],
  },
  {
    type: "interaction.wire-rack", version: 1,
    name: "Snare wires", role: "interaction",
    inputs: [{ name: "motion", type: "audio" }],
    outputs: [{ name: "audio", type: "audio" }],
  },
  {
    type: "observation.three-source", version: 1,
    name: "Observation", role: "observation",
    inputs: [
      { name: "direct", type: "audio" }, { name: "body", type: "audio" },
      { name: "wires", type: "audio" },
    ],
    outputs: [{ name: "audio", type: "audio" }],
  },
];

export const ModuleTypes = new Map(types.map(item => [item.type, item]));

export const RecipeTypes = new Map([
  ["metal.cymbal.v1", { name: "Metallic plate", available: true }],
  ["metal.pair.v1", { name: "Interacting metallic plates", available: false }],
  ["drum.membrane.v1", { name: "Membrane", available: true }],
  ["drum.snare.v1", { name: "Membrane with wires", available: true }],
  ["drum.kick-fm.v1", { name: "Compact FM kick", available: true }],
  ["drum.kick-acoustic.v1", {
    name: "Exciter with resonant body", available: false,
  }],
]);

export function modulePresentation(patch) {
  return new Map(patch.nodes.map(node => {
    const type = ModuleTypes.get(node.type);
    return [node.id, {
      ...type, name: node.name ?? type.name,
      colour: ModuleRoles[type.role].colour,
    }];
  }));
}
