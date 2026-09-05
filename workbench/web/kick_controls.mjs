import {KickModalControls} from "./kick_modal_controls.mjs";

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
  ["kick-thump-controls", [
    "thump_level", "thump_pitch_hz", "thump_pitch_drop_octaves",
    "thump_pitch_fall_seconds", "thump_decay_seconds",
  ]],
  ["kick-resonance-controls", ["resonance_level", "resonance_decay_seconds", "resonance_decay_tilt"]],
  ["kick-tension-controls", ["tension_octaves", "tension_recovery_seconds"]],
  ["kick-contact-controls", [
    "contact_level", "contact_width_seconds", "contact_colour", "contact_noise_level", "contact_noise_decay_seconds",
  ]],
  ["kick-observation-controls", ["equalizer_mode"]],
  ["kick-radiation-controls", ["low_cut_hz", "high_cut_hz", "colour_frequency_hz", "colour_gain_db"]],
  ["kick-multiband-controls", [1,2,3,4].flatMap(n => [`band_${n}_frequency_hz`, `band_${n}_gain_db`])],
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
    this.destroy();
    document.querySelectorAll("[data-kick-controls]").forEach(
      element => element.replaceChildren(),
    );
    for (const [parentId, keys] of Groups)
      for (const key of keys) this.#slider(parentId, key);
    this.#paintEq();
    this.modal = new KickModalControls(this);
  }

  destroy() { this.modal?.destroy(); }
  addSlider(parent, key) { this.#slider(parent, key); }

  refresh(key) {
    const descriptor = this.byKey.get(key);
    const row = document.querySelector(`[data-kick-key="${key}"]`);
    if (!descriptor || !row) return;
    const value = this.state.macros[descriptor.index];
    if (key === "equalizer_mode") {
      row.querySelector("select").value = value;
      this.#paintEq();
      return;
    }
    row.querySelector("input").value = normalized(descriptor, value);
    row.querySelector("output").textContent = valueText(descriptor, value);
  }

  #slider(parentId, key) {
    const descriptor = this.byKey.get(key);
    if (!descriptor) throw new Error(`missing kick parameter ${key}`);
    if (key === "equalizer_mode") { this.#eqChoice(parentId, descriptor); return; }
    const row = document.createElement("label");
    row.className = "slider-row";
    row.dataset.kickKey = key;
    row.dataset.tooltip = `${descriptor.name}; double-click to restore its default.`;
    if (key.endsWith("decay_seconds"))
      row.dataset.tooltip = `${descriptor.name}: time to lose 60 dB of amplitude. Contact decay is the base before implement/spread shaping. Double-click resets.`;
    if (key === "resonance_level")
      row.dataset.tooltip = "Audible resonance only: does not change contact, thump, damping or drive energy. Zero removes the resonance observation. Double-click resets.";
    if (key === "resonance_decay_seconds")
      row.dataset.tooltip = "T60 at 100 Hz. The shared damping slope determines decay at every painted frequency; there are no per-mode decay multipliers.";
    if (key === "resonance_decay_tilt")
      row.dataset.tooltip = "Octaves of T60 lost per frequency octave. +1 halves T60 each octave; 0 is flat. Independent of how many modes exist or their editor slots. Safety bounds: 2 ms to 30 s.";
    if (key === "thump_pitch_fall_seconds")
      row.dataset.tooltip = "Time to reach the settled thump pitch, independent of amplitude decay. Double-click resets.";
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

  #paintEq() {
    const mode = this.state.macros[this.byKey.get("equalizer_mode").index];
    document.getElementById("kick-radiation-controls").hidden = mode === 0;
    for (const key of ["colour_frequency_hz", "colour_gain_db"])
      document.querySelector(`[data-kick-key="${key}"]`).hidden = mode !== 1;
    document.getElementById("kick-multiband-controls").hidden = mode !== 2;
  }

  #eqChoice(parentId, descriptor) {
    const row = document.createElement("label");
    row.dataset.kickKey = descriptor.key;
    row.textContent = "Output EQ ";
    const select = document.createElement("select");
    ["Bypass", "Radiation", "Multiband"].forEach((label, value) => {
      const option = document.createElement("option");
      option.value = value; option.textContent = label; select.append(option);
    });
    select.value = this.state.macros[descriptor.index];
    select.onchange = () => {
      this.state.macros[descriptor.index] = Number(select.value);
      this.#paintEq(); this.onChange(descriptor.key);
    };
    select.ondblclick = e => { e.preventDefault(); select.value = descriptor.defaultValue; select.onchange(); };
    row.append(select); document.getElementById(parentId).append(row);
  }
}
