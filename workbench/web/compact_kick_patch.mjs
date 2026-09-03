import { PatchSchema } from "./percussion_patch.mjs";
import { ModuleTypes } from "./percussion_registry.mjs";

const Nodes = [
  ["kick-primary", "exciter.correlated-fm", 24, 18, "Primary FM"],
  ["kick-secondary", "exciter.correlated-fm", 24, 82, "Secondary FM"],
  ["kick-click", "exciter.noise-burst", 24, 146, "Click"],
  ["kick-mix", "transform.sum3", 250, 82, "Source mix"],
  ["kick-observation", "observation.single-source", 470, 82, "Observation"],
  ["kick-output", "output.mono", 690, 82, "Mono"],
];

const Connections = [
  ["kick-primary.audio", "kick-mix.a"],
  ["kick-secondary.audio", "kick-mix.b"],
  ["kick-click.audio", "kick-mix.c"],
  ["kick-mix.audio", "kick-observation.source"],
  ["kick-observation.audio", "kick-output.audio"],
];

function owner(key) {
  if (key === "model_level_db") return "kick-output";
  if (key.startsWith("secondary_")) return "kick-secondary";
  if (key.startsWith("click_")) return "kick-click";
  if (key === "low_cut_hz" || key === "high_cut_hz")
    return "kick-observation";
  return "kick-primary";
}

function clone(value) { return JSON.parse(JSON.stringify(value)); }

export function createKickPatch(descriptors, values) {
  const patch = {
    schema: PatchSchema,
    id: "factory.kick.compact-01",
    name: "Compact FM kick",
    engineMinimum: 1,
    recipe: "drum.kick-fm.v1",
    nodes: Nodes.map(([id, type, x, y, name]) => ({
      id, type, name, version: ModuleTypes.get(type).version, parameters: {},
      editor: { x, y },
    })),
    connections: Connections.map(([from, to], index) => ({
      id: `kick-route-${index + 1}`, from, to, enabled: true, gain: 1,
      required: index >= 3,
    })),
    outputs: { mono: "kick-output.audio" },
    performanceControls: ["strength", "hardness"],
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
  if (patch.recipe !== "drum.kick-fm.v1" || patch.nodes.length !== Nodes.length ||
      patch.connections.length !== Connections.length ||
      patch.outputs?.mono !== "kick-output.audio") {
    throw new Error("patch is not supported by the compact kick recipe");
  }
  const types = new Map(patch.nodes.map(node => [node.id, node.type]));
  for (const [id, type] of Nodes) {
    if (types.get(id) !== type)
      throw new Error("patch is not supported by the compact kick recipe");
  }
  const routes = new Map(patch.connections.map(item => [
    `${item.from}>${item.to}`, item,
  ]));
  for (const [from, to] of Connections) {
    if (!routes.has(`${from}>${to}`))
      throw new Error("patch is not supported by the compact kick recipe");
  }
  for (const node of patch.nodes) {
    for (const key of Object.keys(node.parameters ?? {})) {
      if (owner(key) !== node.id)
        throw new Error(`invalid kick parameter owner: ${node.id}.${key}`);
    }
  }
  for (const [from, to] of Connections.slice(3)) {
    const route = routes.get(`${from}>${to}`);
    if (route.enabled === false || route.gain !== undefined && route.gain !== 1)
      throw new Error("compact kick output routes are required");
  }
  const audible = Connections.slice(0, 3).some(([from, to]) => {
    const route = routes.get(`${from}>${to}`);
    return route.enabled !== false && route.gain !== 0;
  });
  if (!audible) throw new Error("the compact kick patch has no audible route");
  return patch;
}

export function kickRoutingValues(patch) {
  return Connections.slice(0, 3).map(([from, to]) => {
    const route = patch.connections.find(item =>
      item.from === from && item.to === to);
    return route?.enabled === false ? 0 : route?.gain ?? 1;
  });
}
