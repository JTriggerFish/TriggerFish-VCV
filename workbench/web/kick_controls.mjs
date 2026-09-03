const clamp = (value, minimum, maximum) =>
  Math.max(minimum, Math.min(maximum, value));

function normalized(descriptor, value) {
  if (descriptor.scale === "logarithmic") {
    return Math.log(value / descriptor.minimum) /
      Math.log(descriptor.maximum / descriptor.minimum);
  }
  return (value - descriptor.minimum) /
    (descriptor.maximum - descriptor.minimum);
}

function denormalized(descriptor, position) {
  if (descriptor.scale === "logarithmic") {
    return descriptor.minimum *
      (descriptor.maximum / descriptor.minimum) ** position;
  }
  return descriptor.minimum + position *
    (descriptor.maximum - descriptor.minimum);
}

function valueText(descriptor, value) {
  const digits = descriptor.unit === "Hz" ? (value < 1000 ? 1 : 0) :
    descriptor.unit === "s" ? (value < .1 ? 3 : 2) : 2;
  return `${Number(value).toFixed(digits)}${descriptor.unit ? ` ${descriptor.unit}` : ""}`;
}

const Groups = new Map([
  ["kick-output-controls", ["model_level_db"]],
  ["kick-primary-controls", [
    "fundamental_hz", "pitch_drop_octaves", "pitch_decay_seconds",
    "body_decay_seconds", "fm_depth_hz", "fm_decay_seconds",
    "fm_roughness_hz",
  ]],
  ["kick-secondary-controls", ["secondary_ratio", "secondary_level"]],
  ["kick-click-controls", [
    "click_level", "click_decay_seconds", "click_tilt_db",
  ]],
  ["kick-observation-controls", ["low_cut_hz", "high_cut_hz"]],
]);

export class KickControls {
  constructor({ descriptors, state, onChange, onLevelReset }) {
    this.descriptors = descriptors;
    this.byKey = new Map(descriptors.map(item => [item.key, item]));
    this.state = state;
    this.onChange = onChange;
    this.onLevelReset = onLevelReset;
  }

  build() {
    document.querySelectorAll("[data-kick-controls]").forEach(
      element => element.replaceChildren(),
    );
    for (const [parentId, keys] of Groups)
      for (const key of keys) this.#slider(parentId, key);
  }

  refresh(key) {
    const descriptor = this.byKey.get(key);
    const row = document.querySelector(`[data-kick-key="${key}"]`);
    if (!descriptor || !row) return;
    const value = this.state.macros[descriptor.index];
    row.querySelector("input").value = normalized(descriptor, value);
    row.querySelector("output").textContent = valueText(descriptor, value);
  }

  #slider(parentId, key) {
    const descriptor = this.byKey.get(key);
    if (!descriptor) throw new Error(`missing kick parameter ${key}`);
    const row = document.createElement("label");
    row.className = "slider-row";
    row.dataset.kickKey = key;
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
      const value = this.state.macros[descriptor.index];
      input.value = normalized(descriptor, value);
      output.textContent = valueText(descriptor, value);
    };
    input.oninput = () => {
      this.state.macros[descriptor.index] = clamp(
        denormalized(descriptor, Number(input.value)),
        descriptor.minimum, descriptor.maximum,
      );
      output.textContent = valueText(
        descriptor, this.state.macros[descriptor.index],
      );
      this.onChange(key);
    };
    input.ondblclick = event => {
      event.preventDefault();
      if (key === "model_level_db") this.onLevelReset();
      else {
        this.state.macros[descriptor.index] = descriptor.defaultValue;
        this.onChange(key);
      }
      paint();
    };
    paint();
    row.append(title, input, output);
    document.getElementById(parentId).append(row);
  }
}
