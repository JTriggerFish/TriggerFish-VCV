import {
  ModuleRoles, ModuleTypes, modulePresentation,
} from "./percussion_registry.mjs";

const SvgNamespace = "http://www.w3.org/2000/svg";
const NodeWidth = 124;
const NodeHeight = 48;

function svg(name, attributes = {}) {
  const element = document.createElementNS(SvgNamespace, name);
  for (const [key, value] of Object.entries(attributes)) {
    element.setAttribute(key, value);
  }
  return element;
}

function endpoint(value) {
  const separator = value.lastIndexOf(".");
  return { node: value.slice(0, separator), port: value.slice(separator + 1) };
}

function portY(node, port, direction) {
  const type = ModuleTypes.get(node.type);
  const ports = direction === "input" ? type.inputs : type.outputs;
  const index = Math.max(0, ports.findIndex(item => item.name === port));
  return node.editor.y + NodeHeight * (index + 1) / (ports.length + 1);
}

function routePath(x1, y1, x2, y2, viaY) {
  if (!Number.isFinite(viaY)) {
    const bend = Math.max(26, .45 * (x2 - x1));
    return `M ${x1} ${y1} C ${x1 + bend} ${y1}, ${x2 - bend} ${y2}, ${x2} ${y2}`;
  }
  return `M ${x1} ${y1} ` +
    `C ${x1 + 24} ${y1}, ${x1 + 24} ${viaY}, ${x1 + 48} ${viaY} ` +
    `L ${x2 - 48} ${viaY} ` +
    `C ${x2 - 24} ${viaY}, ${x2 - 24} ${y2}, ${x2} ${y2}`;
}

export class RoutingView {
  constructor({ compact, expanded, dialog, onLayoutChange, onRouteToggle }) {
    this.compact = compact;
    this.expanded = expanded;
    this.dialog = dialog;
    this.onLayoutChange = onLayoutChange;
    this.onRouteToggle = onRouteToggle;
    compact.ondblclick = event => {
      event.preventDefault();
      this.dialog.showModal();
      this.renderExpanded();
    };
  }

  setPatch(patch) {
    this.patch = patch;
    patch.nodes.forEach((node, index) => {
      if (!node.editor) node.editor = {
        x: 24 + 160 * (index % 5), y: 22 + 70 * Math.floor(index / 5),
      };
    });
    this.presentation = modulePresentation(patch);
    this.render(this.compact, false);
    if (this.dialog.open) this.renderExpanded();
    this.decorateSections();
  }

  renderExpanded() {
    this.render(this.expanded, true);
  }

  render(target, editable) {
    target.replaceChildren();
    target.setAttribute("viewBox", "0 0 840 210");
    const nodes = new Map(this.patch.nodes.map(node => [node.id, node]));
    const edges = svg("g", { class: "routing-edges" });
    for (const connection of this.patch.connections) {
      const from = endpoint(connection.from);
      const to = endpoint(connection.to);
      const source = nodes.get(from.node);
      const destination = nodes.get(to.node);
      const x1 = source.editor.x + NodeWidth;
      const y1 = portY(source, from.port, "output");
      const x2 = destination.editor.x;
      const y2 = portY(destination, to.port, "input");
      const path = svg("path", {
        d: routePath(x1, y1, x2, y2, connection.editor?.viaY),
        "data-connection-id": connection.id,
        class: "routing-edge" +
          (connection.enabled === false ? " disabled" : "") +
          (connection.required ? " locked" : "") +
          (editable ? " editable" : ""),
        style: `--route-colour:${this.presentation.get(from.node).colour}`,
      });
      const title = svg("title");
      title.textContent = `${connection.from} → ${connection.to}` +
        (connection.required ? " (required)" : " (click to toggle)");
      path.append(title);
      if (editable && !connection.required) {
        path.onclick = event => {
          event.stopPropagation();
          this.onRouteToggle?.(connection.id);
        };
      }
      edges.append(path);
    }
    target.append(edges);

    for (const node of this.patch.nodes) {
      const type = this.presentation.get(node.id);
      const group = svg("g", {
        class: `routing-node${editable ? " editable" : ""}`,
        transform: `translate(${node.editor.x} ${node.editor.y})`,
        tabindex: "0",
        role: "button",
        "aria-label": `${type.name}, ${ModuleRoles[type.role].name}`,
      });
      const box = svg("rect", {
        width: NodeWidth, height: NodeHeight, rx: 7,
        style: `--module-colour:${type.colour}`,
      });
      const name = svg("text", { x: 12, y: 21, class: "routing-node-name" });
      name.textContent = type.name;
      const role = svg("text", { x: 12, y: 37, class: "routing-node-role" });
      role.textContent = ModuleRoles[type.role].name;
      group.append(box, name, role);
      group.onclick = () => this.select(node.id);
      if (editable) this.bindDrag(group, node);
      target.append(group);
    }
  }

  select(nodeId) {
    const sections = document.querySelectorAll(`[data-module-id="${nodeId}"]`);
    sections.forEach(section => {
      section.classList.remove("module-flash");
      requestAnimationFrame(() => section.classList.add("module-flash"));
    });
    if (!this.dialog.open) sections[0]?.scrollIntoView({ block: "nearest" });
  }

  bindDrag(group, node) {
    let drag;
    group.onpointerdown = event => {
      if (event.button !== 0) return;
      event.preventDefault();
      drag = { id: event.pointerId, x: event.clientX, y: event.clientY,
        nodeX: node.editor.x, nodeY: node.editor.y };
      group.setPointerCapture(event.pointerId);
    };
    group.onpointermove = event => {
      if (!drag || event.pointerId !== drag.id) return;
      const bounds = this.expanded.getBoundingClientRect();
      const scaleX = 840 / bounds.width;
      const scaleY = 210 / bounds.height;
      node.editor.x = Math.max(8, Math.min(708,
        drag.nodeX + (event.clientX - drag.x) * scaleX));
      node.editor.y = Math.max(8, Math.min(154,
        drag.nodeY + (event.clientY - drag.y) * scaleY));
      group.setAttribute("transform", `translate(${node.editor.x} ${node.editor.y})`);
    };
    const finish = event => {
      if (!drag || event.pointerId !== drag.id) return;
      drag = null;
      if (group.hasPointerCapture(event.pointerId))
        group.releasePointerCapture(event.pointerId);
      this.onLayoutChange?.(this.patch);
      this.renderExpanded();
      this.render(this.compact, false);
    };
    group.onpointerup = finish;
    group.onpointercancel = finish;
  }

  decorateSections() {
    document.querySelectorAll("[data-module-id]").forEach(section => {
      const type = this.presentation.get(section.dataset.moduleId);
      section.hidden = !type;
      if (type) section.style.setProperty("--module-colour", type.colour);
    });
  }
}
