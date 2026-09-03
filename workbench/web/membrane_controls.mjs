import {
  membranePresetId, membranePresetName, membranePresetValues,
} from "./membrane_patch.mjs";

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

function valueText(descriptor, value) {
  const digits = descriptor.unit === "Hz" ? (value < 1000 ? 1 : 0) :
    descriptor.unit === "s" ? (value < .1 ? 3 : 2) : 2;
  return `${Number(value).toFixed(digits)}${
    descriptor.unit && descriptor.unit !== "mode" ? ` ${descriptor.unit}` : ""}`;
}

const Groups = new Map([
  ["membrane-output-controls", ["model_level_db"]],
  ["membrane-body-controls", [
    "fundamental_hz", "decay_seconds", "decay_tilt", "inharmonicity",
    "body_brightness",
  ]],
  ["membrane-tension-controls", ["tension_octaves", "tension_decay_seconds"]],
  ["membrane-contact-controls", [
    "contact_level", "contact_duration_seconds", "contact_brightness",
  ]],
  ["membrane-fm-controls", [
    "fm_level", "fm_depth_hz", "fm_decay_seconds", "pitch_drop_octaves",
  ]],
  ["membrane-observation-controls", [
    "direct_level", "body_level", "direct_delay_ms",
  ]],
  ["membrane-radiation-controls", [
    "low_cut_hz", "high_cut_hz", "colour_frequency_hz", "colour_gain_db",
  ]],
  ["membrane-multiband-controls", [
    "band_1_frequency_hz", "band_1_gain_db", "band_2_frequency_hz",
    "band_2_gain_db", "band_3_frequency_hz", "band_3_gain_db",
    "band_4_frequency_hz", "band_4_gain_db",
  ]],
]);

export class MembraneControls {
  constructor({ descriptors, state, onChange, onLevelReset,
                showPresets = true }) {
    this.descriptors = descriptors;
    this.byKey = new Map(descriptors.map(item => [item.key, item]));
    this.state = state;
    this.onChange = onChange;
    this.onLevelReset = onLevelReset;
    this.showPresets = showPresets;
  }

  build() {
    document.querySelectorAll("[data-membrane-controls]").forEach(
      element => element.replaceChildren(),
    );
    if (this.showPresets) this.#presets();
    this.#mode();
    for (const [parentId, keys] of Groups)
      for (const key of keys) this.#slider(parentId, key);
    this.#paintMode();
  }

  refresh(key) {
    const descriptor = this.byKey.get(key);
    const row = document.querySelector(`[data-membrane-key="${key}"]`);
    if (!descriptor || !row) return;
    const value = this.state.macros[descriptor.index];
    row.querySelector("input").value = normalized(descriptor, value);
    row.querySelector("output").textContent = valueText(descriptor, value);
  }

  #presets() {
    const parent = document.getElementById("membrane-preset-controls");
    for (const key of ["tom", "acousticKick"]) {
      const button = document.createElement("button");
      button.type = "button";
      button.textContent = membranePresetName(key);
      button.dataset.tooltip = `Load the ${button.textContent} starting point; routing remains editable.`;
      button.onclick = () => {
        const values = membranePresetValues(key, this.descriptors);
        this.state.macros.splice(0, this.state.macros.length, ...values);
        this.state.patch.id = membranePresetId(key);
        this.state.patch.name = membranePresetName(key);
        this.build();
        this.onChange("preset");
      };
      parent.append(button);
    }
  }

  #mode() {
    const descriptor = this.byKey.get("equalizer_mode");
    const parent = document.getElementById("membrane-eq-mode-controls");
    const label = document.createElement("label");
    label.textContent = "Output EQ ";
    label.dataset.tooltip = "Bypass, the compact radiation filter, or four independent parametric bands.";
    const select = document.createElement("select");
    for (const [value, name] of [[0, "Bypass"], [1, "Radiation"], [2, "Multiband"]]) {
      const option = document.createElement("option");
      option.value = value;
      option.textContent = name;
      select.append(option);
    }
    select.value = Math.round(this.state.macros[descriptor.index]);
    select.onchange = () => {
      this.state.macros[descriptor.index] = Number(select.value);
      this.#paintMode();
      this.onChange("equalizer_mode");
    };
    select.ondblclick = event => {
      event.preventDefault();
      select.value = descriptor.defaultValue;
      select.dispatchEvent(new Event("change"));
    };
    label.append(select);
    parent.append(label);
  }

  #paintMode() {
    const descriptor = this.byKey.get("equalizer_mode");
    const mode = Math.round(this.state.macros[descriptor.index]);
    document.getElementById("membrane-radiation-panel").hidden = mode !== 1;
    document.getElementById("membrane-multiband-panel").hidden = mode !== 2;
  }

  #slider(parentId, key) {
    const descriptor = this.byKey.get(key);
    if (!descriptor) throw new Error(`missing membrane parameter ${key}`);
    const row = document.createElement("label");
    row.className = "slider-row";
    row.dataset.membraneKey = key;
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
