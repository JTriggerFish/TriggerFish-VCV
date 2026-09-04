import { ModuleTypes, RecipeTypes } from "./percussion_registry.mjs";

export const PatchSchema = "triggerfish.percussion.patch/v1";
export const PatchEngineVersion = 1;
const PerformanceControls = new Set([
  "strength", "location", "hardness", "implement", "contactSpread",
  "mute", "pedal",
]);

export function endpoint(value) {
  const separator = value.lastIndexOf(".");
  if (separator <= 0 || separator + 1 >= value.length) return null;
  return { node: value.slice(0, separator), port: value.slice(separator + 1) };
}

function port(type, direction, name) {
  return type[direction].find(item => item.name === name);
}

function validateNode(node, nodes, descriptorByKey) {
  const type = ModuleTypes.get(node?.type);
  if (!type || type.version !== node.version || typeof node.id !== "string" ||
      !/^[a-z][a-z0-9_-]*$/.test(node.id) || nodes.has(node.id)) {
    throw new Error(`invalid percussion module: ${node?.id ?? "unknown"}`);
  }
  for (const [key, parameter] of Object.entries(node.parameters ?? {})) {
    const descriptor = descriptorByKey.get(key);
    if (!descriptor || !Number.isFinite(parameter) ||
        descriptor.scale === "choice" && !Number.isInteger(parameter) ||
        parameter < descriptor.minimum || parameter > descriptor.maximum) {
      throw new Error(`invalid module parameter: ${node.id}.${key}`);
    }
  }
  if (node.editor && (!Number.isFinite(node.editor.x) ||
      !Number.isFinite(node.editor.y))) {
    throw new Error(`invalid module layout: ${node.id}`);
  }
  nodes.set(node.id, node);
}

function connect(connection, nodes, connectionIds, adjacency, indegree) {
  const from = endpoint(connection?.from ?? "");
  const to = endpoint(connection?.to ?? "");
  const fromNode = from && nodes.get(from.node);
  const toNode = to && nodes.get(to.node);
  const fromType = fromNode && ModuleTypes.get(fromNode.type);
  const toType = toNode && ModuleTypes.get(toNode.type);
  const output = fromType && port(fromType, "outputs", from.port);
  const input = toType && port(toType, "inputs", to.port);
  const invalidEnabled = typeof connection.enabled !== "boolean";
  const invalidGain = Object.hasOwn(connection, "gain");
  if (!output || !input || output.type !== input.type || invalidEnabled ||
      invalidGain || typeof connection.id !== "string" ||
      connectionIds.has(connection.id)) {
    throw new Error(`invalid percussion connection: ${connection?.id ?? "unknown"}`);
  }
  connectionIds.add(connection.id);
  if (connection.enabled === false) return;
  adjacency.get(from.node).push(to.node);
  indegree.set(to.node, indegree.get(to.node) + 1);
}

function validateAcyclic(nodes, adjacency, indegree) {
  const ready = [...nodes.keys()].filter(id => indegree.get(id) === 0);
  let visited = 0;
  while (ready.length) {
    const id = ready.pop();
    ++visited;
    for (const target of adjacency.get(id)) {
      indegree.set(target, indegree.get(target) - 1);
      if (indegree.get(target) === 0) ready.push(target);
    }
  }
  if (visited !== nodes.size)
    throw new Error("percussion patch graph contains a cycle");
}

function validateOutputs(outputs, nodes) {
  if (!outputs.length) throw new Error("percussion patch has no output");
  for (const target of outputs) {
    const parsed = endpoint(target);
    const node = parsed && nodes.get(parsed.node);
    if (!node || !port(ModuleTypes.get(node.type), "outputs", parsed.port))
      throw new Error(`invalid percussion output: ${target}`);
  }
}

function validatePerformanceControls(controls) {
  if (!Array.isArray(controls) || new Set(controls).size !== controls.length ||
      controls.some(control => !PerformanceControls.has(control))) {
    throw new Error("invalid percussion performance controls");
  }
}

export function validatePatch(value, descriptors) {
  const recipe = RecipeTypes.get(value?.recipe);
  if (value?.schema !== PatchSchema || typeof value.id !== "string" ||
      typeof value.name !== "string" || !Number.isInteger(value.engineMinimum) ||
      value.engineMinimum < 1 || value.engineMinimum > PatchEngineVersion ||
      !recipe || !recipe.available || !Array.isArray(value.nodes) ||
      !Array.isArray(value.connections) || value.nodes.length < 1 ||
      value.nodes.length > 32 || value.connections.length > 64) {
    throw new Error("unsupported or incomplete percussion patch");
  }

  const descriptorByKey = new Map(descriptors.map(item => [item.key, item]));
  if (descriptorByKey.size !== descriptors.length)
    throw new Error("percussion parameter descriptors are not unique");
  const nodes = new Map();
  for (const node of value.nodes) validateNode(node, nodes, descriptorByKey);
  const parameterKeys = value.nodes.flatMap(node =>
    Object.keys(node.parameters ?? {}));
  const present = new Set(parameterKeys);
  if (present.size !== parameterKeys.length || present.size !== descriptors.length ||
      descriptors.some(descriptor => !present.has(descriptor.key))) {
    throw new Error("percussion patch parameter set is incomplete");
  }
  const adjacency = new Map([...nodes.keys()].map(id => [id, []]));
  const indegree = new Map([...nodes.keys()].map(id => [id, 0]));
  const connectionIds = new Set();
  for (const connection of value.connections)
    connect(connection, nodes, connectionIds, adjacency, indegree);
  validateAcyclic(nodes, adjacency, indegree);
  validateOutputs(Object.values(value.outputs ?? {}), nodes);
  validatePerformanceControls(value.performanceControls);
  return value;
}
