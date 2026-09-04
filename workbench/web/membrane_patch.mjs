import { PatchSchema } from "./percussion_patch.mjs";
import { ModuleTypes } from "./percussion_registry.mjs";

const Nodes = [
  ["membrane-contact", "exciter.contact", 20, 22, "Contact"],
  ["membrane-fm", "exciter.correlated-fm", 20, 132, "FM supplement"],
  ["membrane-direct-mix", "transform.sum2", 230, 20, "Direct mix"],
  ["membrane-body-mix", "transform.sum2", 230, 132, "Body drive"],
  ["membrane-tension", "interaction.strike-energy", 442, 178, "Strike / tension"],
  ["membrane-body", "body.membrane-modal", 442, 90, "Membrane"],
  ["membrane-observation", "observation.dual-source", 660, 62, "Observation"],
  ["membrane-eq", "observation.equalizer", 856, 62, "Output EQ"],
  ["membrane-output", "output.mono", 1040, 62, "Mono"],
];

const Connections = [
  ["membrane-contact.direct", "membrane-direct-mix.a", "contact-direct", false],
  ["membrane-contact.body", "membrane-body-mix.a", "contact-body", false],
  ["membrane-fm.audio", "membrane-direct-mix.b", "fm-direct", false],
  ["membrane-fm.audio", "membrane-body-mix.b", "fm-body", false],
  ["membrane-body.audio", "membrane-observation.body", "body-observation", false],
  ["membrane-contact.event", "membrane-tension.strike", "strike-energy", true],
  ["membrane-body-mix.audio", "membrane-body.drive", "body-drive", true],
  ["membrane-tension.tension", "membrane-body.tension", "tension", true],
  ["membrane-direct-mix.audio", "membrane-observation.direct", "direct-observation", true],
  ["membrane-observation.audio", "membrane-eq.audio", "observation-eq", true],
  ["membrane-eq.audio", "membrane-output.audio", "eq-output", true],
];

const RoutedEndpoints = [
  ["membrane-contact.direct", "membrane-direct-mix.a"],
  ["membrane-contact.body", "membrane-body-mix.a"],
  ["membrane-fm.audio", "membrane-direct-mix.b"],
  ["membrane-fm.audio", "membrane-body-mix.b"],
  ["membrane-body.audio", "membrane-observation.body"],
];

function owner(key) {
  if (key === "model_level_db") return "membrane-output";
  if (key === "contact_direct_level" || key === "fm_direct_level")
    return "membrane-direct-mix";
  if (key === "contact_body_level" || key === "fm_body_level")
    return "membrane-body-mix";
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

function clone(value) { return JSON.parse(JSON.stringify(value)); }

const Presets = Object.freeze({
  tom: {
    id: "factory.membrane.tom-01",
    name: "Open tom",
    values: {
      fundamental_hz: 105, decay_seconds: 1.15, decay_tilt: .55,
      inharmonicity: .35, body_brightness: .55, tension_octaves: .11,
      tension_decay_seconds: .13,
      contact_direct_level: .245, contact_body_level: .7,
      contact_duration_seconds: .004, contact_brightness: .58,
      fm_direct_level: .0144, fm_body_level: .081,
      fm_depth_hz: 260, fm_decay_seconds: .07,
      pitch_drop_octaves: .28, direct_level: .3, body_level: 1,
      low_cut_hz: 24, high_cut_hz: 18000,
    },
  },
  acousticKick: {
    id: "factory.membrane.acoustic-kick-01",
    name: "Acoustic kick",
    values: {
      model_level_db: -12,
      fundamental_hz: 35, decay_seconds: .25, decay_tilt: 0,
      inharmonicity: .18, body_brightness: 1, tension_octaves: .1,
      tension_decay_seconds: .05,
      contact_direct_level: 1.64, contact_body_level: .205,
      contact_duration_seconds: .01, contact_brightness: 0,
      fm_direct_level: .04, fm_body_level: .05625,
      fm_depth_hz: 520, fm_decay_seconds: .07,
      pitch_drop_octaves: 1, direct_level: 2, body_level: 1.15,
      low_cut_hz: 16, high_cut_hz: 3500,
      colour_frequency_hz: 600, colour_gain_db: 10,
    },
  },
});

export function createMembranePatch(descriptors, values) {
  return createTomPatch(descriptors, values);
}

function createPatch(descriptors, values, preset) {
  const patch = {
    schema: PatchSchema,
    id: preset.id,
    name: preset.name,
    engineMinimum: 1,
    recipe: "drum.membrane.v1",
    nodes: Nodes.map(([id, type, x, y, name]) => ({
      id, type, name, version: ModuleTypes.get(type).version, parameters: {},
      editor: { x, y },
    })),
    connections: Connections.map(([from, to, suffix, required]) => ({
      id: `membrane-route-${suffix}`, from, to, enabled: true,
      required,
    })),
    outputs: { mono: "membrane-output.audio" },
    performanceControls: [
      "strength", "location", "hardness", "implement", "contactSpread",
    ],
  };
  return patchWithMembraneValues(patch, descriptors, values);
}

export function createTomPatch(descriptors, values) {
  return createPatch(descriptors, values, Presets.tom);
}

export function createAcousticKickPatch(descriptors) {
  const values = membranePresetValues("acousticKick", descriptors);
  return createPatch(descriptors, values, Presets.acousticKick);
}

export function patchWithMembraneValues(patch, descriptors, values) {
  const result = clone(patch);
  const nodes = new Map(result.nodes.map(node => [node.id, node]));
  descriptors.forEach((descriptor, index) => {
    nodes.get(owner(descriptor.key)).parameters[descriptor.key] = values[index];
  });
  return result;
}

export function membraneValuesFromPatch(patch, descriptors) {
  const values = descriptors.map(item => item.defaultValue);
  const indices = new Map(descriptors.map((item, index) => [item.key, index]));
  for (const node of patch.nodes) {
    for (const [key, value] of Object.entries(node.parameters ?? {})) {
      const index = indices.get(key);
      if (index !== undefined) values[index] = value;
    }
  }
  return values;
}

export function membranePresetValues(name, descriptors) {
  const preset = Presets[name];
  if (!preset) throw new Error(`unknown membrane preset: ${name}`);
  return descriptors.map(descriptor =>
    preset.values[descriptor.key] ?? descriptor.defaultValue);
}

export function membranePresetName(name) { return Presets[name]?.name; }
export function membranePresetId(name) { return Presets[name]?.id; }

export function validateMembranePatch(patch) {
  if (patch.recipe !== "drum.membrane.v1" ||
      patch.nodes.length !== Nodes.length ||
      patch.connections.length !== Connections.length ||
      patch.outputs?.mono !== "membrane-output.audio")
    throw new Error("patch is not supported by the membrane recipe");
  const types = new Map(patch.nodes.map(node => [node.id, node.type]));
  for (const [id, type] of Nodes)
    if (types.get(id) !== type)
      throw new Error("patch is not supported by the membrane recipe");
  const routes = new Map(patch.connections.map(item => [
    `${item.from}>${item.to}`, item,
  ]));
  for (const [from, to, , required] of Connections) {
    const route = routes.get(`${from}>${to}`);
    if (!route) throw new Error("patch is not supported by the membrane recipe");
    if (required && route.enabled === false)
      throw new Error("membrane structural routes are required");
  }
  for (const node of patch.nodes)
    for (const key of Object.keys(node.parameters ?? {}))
      if (owner(key) !== node.id)
        throw new Error(`invalid membrane parameter owner: ${node.id}.${key}`);
  const active = ([from, to]) => {
    const route = routes.get(`${from}>${to}`);
    return route?.enabled !== false;
  };
  const direct = active(RoutedEndpoints[0]) || active(RoutedEndpoints[2]);
  const body = active(RoutedEndpoints[4]) &&
    (active(RoutedEndpoints[1]) || active(RoutedEndpoints[3]));
  if (!direct && !body)
    throw new Error("the membrane patch has no audible route");
  return patch;
}

export function membraneRoutingValues(patch) {
  return RoutedEndpoints.map(([from, to]) => {
    const route = patch.connections.find(item =>
      item.from === from && item.to === to);
    return route?.enabled !== false;
  });
}
