import { PatchSchema, endpoint } from "./percussion_patch.mjs";
import { ModuleTypes } from "./percussion_registry.mjs";

const Nodes = [
  ["contact", "exciter.contact", 24, 92],
  ["bloom", "transform.dispersion-loop", 184, 22],
  ["body", "body.stochastic-modal-field", 350, 92],
  ["observation", "observation.dual-source", 520, 92],
  ["output", "output.mono", 690, 92],
];

const Connections = [
  ["contact.body", "bloom.drive"],
  ["contact.body", "body.primary"],
  ["bloom.return", "body.secondary"],
  ["contact.direct", "observation.direct"],
  ["body.audio", "observation.body"],
  ["observation.audio", "output.audio"],
];

const TypeByNode = Object.fromEntries(
  Nodes.map(([id, type]) => [id, type]),
);

function moduleForParameter(key) {
  if (key === "model_level_db") return "output";
  if (key.startsWith("direct_radiation_") || key.startsWith("direct_low_") ||
      key.startsWith("direct_high_") || key.startsWith("direct_colour_") ||
      key.startsWith("dense_radiation_") || key.startsWith("dense_low_") ||
      key.startsWith("dense_high_") || key.startsWith("dense_colour_")) {
    return "observation";
  }
  if (key === "direct_gain" || key.startsWith("impact_")) return "contact";
  if (key.startsWith("bloom_")) return "bloom";
  if (key === "field_gain" || key === "body_brightness" ||
      key.startsWith("field_") || key.startsWith("body_decay_") ||
      key.startsWith("resolved_")) return "body";
  return null;
}

function clone(value) {
  return JSON.parse(JSON.stringify(value));
}

export function createCrashPatch(descriptors, values) {
  const nodes = Nodes.map(([id, type, x, y]) => ({
    id, type, version: ModuleTypes.get(type).version, parameters: {},
    editor: { x, y },
  }));
  const patch = {
    schema: PatchSchema,
    id: "factory.crash.experimental-01",
    name: "Experimental crash",
    engineMinimum: 1,
    recipe: "metal.cymbal.v1",
    nodes,
    connections: Connections.map(([from, to], index) => ({
      id: `route-${index + 1}`, from, to, enabled: true, gain: 1,
      required: index + 1 === Connections.length,
      ...(from === "contact.direct" ? { editor: { viaY: 184 } } : {}),
    })),
    outputs: { mono: "output.audio" },
    performanceControls: [
      "strength", "location", "hardness", "implement", "contactSpread",
    ],
  };
  return patchWithMacroValues(patch, descriptors, values);
}

export function patchWithMacroValues(patch, descriptors, values) {
  const result = clone(patch);
  const nodes = new Map(result.nodes.map(node => [node.id, node]));
  descriptors.forEach((descriptor, index) => {
    const owner = moduleForParameter(descriptor.key);
    if (owner && nodes.has(owner))
      nodes.get(owner).parameters[descriptor.key] = values[index];
  });
  return result;
}

export function macroValuesFromPatch(patch, descriptors) {
  const values = descriptors.map(item => item.defaultValue);
  const byKey = new Map(descriptors.map((item, index) => [item.key, index]));
  for (const node of patch.nodes) {
    for (const [key, value] of Object.entries(node.parameters ?? {})) {
      const index = byKey.get(key);
      if (index !== undefined) values[index] = value;
    }
  }
  return values;
}

function validateParameterOwnership(patch) {
  for (const node of patch.nodes) {
    for (const key of Object.keys(node.parameters ?? {})) {
      if (moduleForParameter(key) !== node.id)
        throw new Error(`invalid module parameter owner: ${node.id}.${key}`);
    }
  }
}

function validateAudibleRoute(patch) {
  const reachable = new Set(["contact"]);
  for (const node of ["bloom", "body", "observation", "output"]) {
    const connected = patch.connections.some(connection => {
      const from = endpoint(connection.from);
      const to = endpoint(connection.to);
      return to.node === node && reachable.has(from.node) &&
        connection.enabled !== false && connection.gain !== 0;
    });
    if (connected) reachable.add(node);
  }
  if (!reachable.has("output"))
    throw new Error("the compiled crash patch has no audible route");
}

export function validateCrashAdapterPatch(patch) {
  if (patch.recipe !== "metal.cymbal.v1" || patch.nodes.length !== Nodes.length ||
      patch.connections.length !== Connections.length ||
      patch.outputs?.mono !== "output.audio") {
    throw new Error("patch is not supported by the metallic-plate recipe");
  }
  const types = new Map(patch.nodes.map(node => [node.id, node.type]));
  for (const [id, type] of Nodes.map(([id, type]) => [id, type])) {
    if (types.get(id) !== type)
      throw new Error("patch is not supported by the metallic-plate recipe");
  }
  const routes = new Set(patch.connections.map(item => `${item.from}>${item.to}`));
  for (const [from, to] of Connections) {
    if (!routes.has(`${from}>${to}`))
      throw new Error("patch is not supported by the metallic-plate recipe");
  }
  const output = patch.connections.find(item =>
    item.from === "observation.audio" && item.to === "output.audio");
  if (output.enabled === false || output.gain !== undefined && output.gain !== 1)
    throw new Error("the compiled crash output route is required");
  validateParameterOwnership(patch);
  validateAudibleRoute(patch);
  return patch;
}

export function routingValuesFromPatch(patch) {
  return Connections.slice(0, -1).map(([from, to]) => {
    const connection = patch.connections.find(item =>
      item.from === from && item.to === to);
    return connection?.enabled === false ? 0 : connection?.gain ?? 1;
  });
}
