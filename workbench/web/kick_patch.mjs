import { PatchSchema } from "./percussion_patch.mjs";
import { ModuleTypes } from "./percussion_registry.mjs";

const Nodes = [
  ["kick-contact", "exciter.contact", 24, 18, "Contact"],
  ["kick-thump", "exciter.thump", 24, 185, "Thump"],
  ["kick-resonance", "body.membrane-modal", 270, 18, "Resonance"],
  ["kick-tension", "interaction.strike-energy", 24, 100, "Strike / tension"],
  ["kick-mix", "transform.sum3", 480, 110, "Observation mix"],
  ["kick-observation", "observation.equalizer", 690, 110, "Output EQ"],
  ["kick-output", "output.mono", 900, 110, "Mono"],
];

const Connections = [
  ["kick-contact.direct", "kick-mix.a"],
  ["kick-thump.audio", "kick-mix.b"],
  ["kick-resonance.audio", "kick-mix.c"],
  ["kick-contact.body", "kick-resonance.drive"],
  ["kick-contact.event", "kick-tension.strike"],
  ["kick-tension.tension", "kick-resonance.tension"],
  ["kick-mix.audio", "kick-observation.audio"],
  ["kick-observation.audio", "kick-output.audio"],
];

function owner(key) {
  if (key === "model_level_db") return "kick-output";
  if (key.startsWith("contact_")) return "kick-contact";
  if (key.startsWith("thump_")) return "kick-thump";
  if (key.startsWith("resonance_")) return "kick-resonance";
  if (key.startsWith("tension_")) return "kick-tension";
  return "kick-observation";
}

function clone(value) { return JSON.parse(JSON.stringify(value)); }

export function createKickPatch(descriptors, values) {
  const patch = {
    schema: PatchSchema,
    id: "factory.kick.standard-01",
    name: "Kick",
    engineMinimum: 1,
    recipe: "drum.kick.v1",
    nodes: Nodes.map(([id, type, x, y, name]) => ({
      id, type, name, version: ModuleTypes.get(type).version, parameters: {},
      editor: { x, y },
    })),
    connections: Connections.map(([from, to], index) => ({
      id: `kick-route-${index + 1}`, from, to, enabled: true,
      required: index >= 3,
    })),
    outputs: { mono: "kick-output.audio" },
    performanceControls: ["strength", "hardness", "implement", "location", "contactSpread"],
  };
  return patchWithKickValues(patch, descriptors, values);
}

export function patchWithKickValues(patch, descriptors, values) {
  const result = clone(patch);
  const nodes = new Map(result.nodes.map(node => [node.id, node]));
  descriptors.forEach((descriptor, index) => {
    nodes.get(owner(descriptor.key)).parameters[descriptor.key] = values[index];
  });
  return result;
}

export function kickValuesFromPatch(patch, descriptors) {
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

export function validateKickPatch(patch) {
  if (patch.recipe !== "drum.kick.v1" || patch.nodes.length !== Nodes.length ||
      patch.connections.length !== Connections.length ||
      patch.outputs?.mono !== "kick-output.audio") {
    throw new Error("patch is not supported by the kick recipe");
  }
  const types = new Map(patch.nodes.map(node => [node.id, node.type]));
  for (const [id, type] of Nodes) {
    if (types.get(id) !== type)
      throw new Error("patch is not supported by the kick recipe");
  }
  const routes = new Map(patch.connections.map(item => [
    `${item.from}>${item.to}`, item,
  ]));
  for (const [from, to] of Connections) {
    if (!routes.has(`${from}>${to}`))
      throw new Error("patch is not supported by the kick recipe");
  }
  for (const node of patch.nodes) {
    for (const key of Object.keys(node.parameters ?? {})) {
      if (owner(key) !== node.id)
        throw new Error(`invalid kick parameter owner: ${node.id}.${key}`);
    }
  }
  for (const [from, to] of Connections.slice(3)) {
    const route = routes.get(`${from}>${to}`);
    if (route.enabled === false)
      throw new Error("kick output routes are required");
  }
  const audible = Connections.slice(0, 3).some(([from, to]) => {
    const route = routes.get(`${from}>${to}`);
    return route.enabled !== false;
  });
  if (!audible) throw new Error("the kick patch has no audible route");
  return patch;
}

export function kickRoutingValues(patch) {
  return Connections.slice(0, 3).map(([from, to]) => {
    const route = patch.connections.find(item =>
      item.from === from && item.to === to);
    return route?.enabled !== false;
  });
}
