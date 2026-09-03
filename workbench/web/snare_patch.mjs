import { PatchSchema } from "./percussion_patch.mjs";
import { ModuleTypes } from "./percussion_registry.mjs";

const Nodes = [
  ["membrane-contact", "exciter.contact", 20, 18, "Contact"],
  ["membrane-fm", "exciter.correlated-fm", 20, 142, "FM supplement"],
  ["membrane-direct-mix", "transform.sum2", 220, 18, "Direct mix"],
  ["membrane-body-mix", "transform.sum2", 220, 142, "Body drive"],
  ["membrane-tension", "interaction.strike-energy", 420, 202, "Strike / tension"],
  ["membrane-body", "body.membrane-modal", 420, 108, "Membrane"],
  ["snare-wires", "interaction.wire-rack", 620, 174, "Snare wires"],
  ["membrane-observation", "observation.three-source", 800, 74, "Observation"],
  ["membrane-eq", "observation.equalizer", 994, 74, "Output EQ"],
  ["membrane-output", "output.mono", 1170, 74, "Mono"],
];

const Connections = [
  ["membrane-contact.direct", "membrane-direct-mix.a", false, .35],
  ["membrane-contact.body", "membrane-body-mix.a", false, 1],
  ["membrane-fm.audio", "membrane-direct-mix.b", false, .05],
  ["membrane-fm.audio", "membrane-body-mix.b", false, .32],
  ["membrane-body.audio", "membrane-observation.body", false, 1],
  ["membrane-body.audio", "snare-wires.motion", false, 1],
  ["snare-wires.audio", "membrane-observation.wires", false, 1],
  ["membrane-contact.event", "membrane-tension.strike", true, 1],
  ["membrane-body-mix.audio", "membrane-body.drive", true, 1],
  ["membrane-tension.tension", "membrane-body.tension", true, 1],
  ["membrane-direct-mix.audio", "membrane-observation.direct", true, 1],
  ["membrane-observation.audio", "membrane-eq.audio", true, 1],
  ["membrane-eq.audio", "membrane-output.audio", true, 1],
];

const RoutedEndpoints = Connections.slice(0, 7).map(item => item.slice(0, 2));
const clone = value => JSON.parse(JSON.stringify(value));

function owner(key) {
  if (key.startsWith("wire_")) return "snare-wires";
  if (key.startsWith("ring_")) return "membrane-body";
  if (key === "model_level_db") return "membrane-output";
  if (key.startsWith("contact_")) return "membrane-contact";
  if (key.startsWith("fm_") || key === "pitch_drop_octaves")
    return "membrane-fm";
  if (key.startsWith("tension_")) return "membrane-tension";
  if (["fundamental_hz", "decay_seconds", "decay_tilt", "inharmonicity",
    "body_brightness"].includes(key)) return "membrane-body";
  if (["direct_level", "body_level", "direct_delay_ms"].includes(key))
    return "membrane-observation";
  return "membrane-eq";
}

export function createSnarePatch(descriptors, values) {
  const patch = {
    schema: PatchSchema, id: "factory.snare.acoustic-01",
    name: "Acoustic snare", engineMinimum: 1, recipe: "drum.snare.v1",
    nodes: Nodes.map(([id, type, x, y, name]) => ({
      id, type, name, version: ModuleTypes.get(type).version,
      parameters: {}, editor: { x, y },
    })),
    connections: Connections.map(([from, to, required, gain], index) => ({
      id: `snare-route-${index}`, from, to, enabled: true, gain, required,
    })),
    outputs: { mono: "membrane-output.audio" },
    performanceControls: [
      "strength", "location", "hardness", "implement", "contactSpread",
    ],
  };
  return patchWithSnareValues(patch, descriptors, values);
}

export function patchWithSnareValues(patch, descriptors, values) {
  const result = clone(patch);
  const nodes = new Map(result.nodes.map(node => [node.id, node]));
  descriptors.forEach((descriptor, index) => {
    nodes.get(owner(descriptor.key)).parameters[descriptor.key] = values[index];
  });
  return result;
}

export function snareValuesFromPatch(patch, descriptors) {
  const values = descriptors.map(item => item.defaultValue);
  const indices = new Map(descriptors.map((item, index) => [item.key, index]));
  for (const node of patch.nodes)
    for (const [key, value] of Object.entries(node.parameters ?? {})) {
      const index = indices.get(key);
      if (index !== undefined) values[index] = value;
    }
  return values;
}

export function snareRoutingValues(patch) {
  return RoutedEndpoints.map(([from, to]) => {
    const route = patch.connections.find(item =>
      item.from === from && item.to === to);
    return route?.enabled === false ? 0 : route?.gain ?? 1;
  });
}

export function validateSnarePatch(patch) {
  if (patch.recipe !== "drum.snare.v1" ||
      patch.nodes.length !== Nodes.length ||
      patch.connections.length !== Connections.length ||
      patch.outputs?.mono !== "membrane-output.audio")
    throw new Error("patch is not supported by the snare recipe");
  const types = new Map(patch.nodes.map(node => [node.id, node.type]));
  for (const [id, type] of Nodes)
    if (types.get(id) !== type)
      throw new Error("patch is not supported by the snare recipe");
  const routes = new Map(patch.connections.map(item => [
    `${item.from}>${item.to}`, item,
  ]));
  for (const [from, to, required] of Connections) {
    const route = routes.get(`${from}>${to}`);
    if (!route) throw new Error("patch is not supported by the snare recipe");
    if (required && (route.enabled === false ||
        route.gain !== undefined && route.gain !== 1))
      throw new Error("snare structural routes are required");
  }
  for (const node of patch.nodes)
    for (const key of Object.keys(node.parameters ?? {}))
      if (owner(key) !== node.id)
        throw new Error(`invalid snare parameter owner: ${node.id}.${key}`);
  const active = index => {
    const [from, to] = RoutedEndpoints[index];
    const route = routes.get(`${from}>${to}`);
    return route?.enabled !== false && route?.gain !== 0;
  };
  const direct = active(0) || active(2);
  const membraneDrive = active(1) || active(3);
  const body = membraneDrive && active(4);
  const wires = membraneDrive && active(5) && active(6);
  if (!direct && !body && !wires)
    throw new Error("the snare patch has no audible route");
  return patch;
}
