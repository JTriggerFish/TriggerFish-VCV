import { MembraneControls } from "./membrane_controls.mjs";

const clamp = (value, minimum, maximum) =>
  Math.max(minimum, Math.min(maximum, value));

function normalized(descriptor, value) {
  if (descriptor.scale === "logarithmic")
    return Math.log(value / descriptor.minimum) /
      Math.log(descriptor.maximum / descriptor.minimum);
  return (value - descriptor.minimum) /
    (descriptor.maximum - descriptor.minimum);
}

function denormalized(descriptor, position) {
  if (descriptor.scale === "logarithmic")
    return descriptor.minimum *
      (descriptor.maximum / descriptor.minimum) ** position;
  return descriptor.minimum + position *
    (descriptor.maximum - descriptor.minimum);
}

function text(descriptor, value) {
  const digits = descriptor.unit === "Hz" ? (value < 1000 ? 1 : 0) :
    descriptor.unit === "s" ? (value < .1 ? 3 : 2) : 2;
  return `${Number(value).toFixed(digits)}${
    descriptor.unit ? ` ${descriptor.unit}` : ""}`;
}

const Groups = new Map([
  ["snare-ring-controls", [
    "ring_frequency_hz", "ring_decay_seconds", "ring_level",
  ]],
  ["snare-wire-response-controls", [
    "wire_level", "wire_sensitivity", "wire_threshold",
    "wire_motion_highpass_hz", "wire_attack_seconds", "wire_release_seconds",
  ]],
  ["snare-wire-spectrum-controls", [
    "wire_minimum_hz", "wire_maximum_hz", "wire_decay_seconds",
    "wire_decay_tilt", "wire_density", "wire_brightness",
    "wire_noise_mix", "wire_modal_mix",
  ]],
]);

export class SnareControls extends MembraneControls {
  constructor(options) {
    super({ ...options, showPresets: false });
    this.snareDescriptors = options.descriptors;
    this.snareState = options.state;
    this.snareOnChange = options.onChange;
  }

  build() {
    super.build();
    document.querySelectorAll("[data-snare-controls]").forEach(
      element => element.replaceChildren(),
    );
    this.#preset();
    for (const [parent, keys] of Groups)
      for (const key of keys) this.#slider(parent, key);
  }

  #preset() {
    const button = document.createElement("button");
    button.type = "button";
    button.textContent = "Acoustic snare";
    button.dataset.tooltip = "Restore the compiled acoustic-snare starting point.";
    button.onclick = () => {
      this.snareState.macros.splice(0, this.snareState.macros.length,
        ...this.snareDescriptors.map(item => item.defaultValue));
      this.build();
      this.snareOnChange("preset");
    };
    document.getElementById("membrane-preset-controls").append(button);
  }

  #slider(parentId, key) {
    const descriptor = this.snareDescriptors.find(item => item.key === key);
    if (!descriptor) throw new Error(`missing snare parameter ${key}`);
    const row = document.createElement("label");
    row.className = "slider-row";
    row.dataset.tooltip = `${descriptor.name}; double-click to restore its default.`;
    const title = document.createElement("span");
    title.textContent = descriptor.name;
    const input = document.createElement("input");
    input.type = "range";
    input.min = 0;
    input.max = 1;
    input.step = 1 / 500;
    const output = document.createElement("output");
    const paint = () => {
      const value = this.snareState.macros[descriptor.index];
      input.value = normalized(descriptor, value);
      output.textContent = text(descriptor, value);
    };
    input.oninput = () => {
      this.snareState.macros[descriptor.index] = clamp(
        denormalized(descriptor, Number(input.value)),
        descriptor.minimum, descriptor.maximum,
      );
      output.textContent = text(
        descriptor, this.snareState.macros[descriptor.index]);
      this.snareOnChange(key);
    };
    input.ondblclick = event => {
      event.preventDefault();
      this.snareState.macros[descriptor.index] = descriptor.defaultValue;
      paint();
      this.snareOnChange(key);
    };
    paint();
    row.append(title, input, output);
    document.getElementById(parentId).append(row);
  }
}
