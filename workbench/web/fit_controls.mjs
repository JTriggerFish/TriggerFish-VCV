import { PointEditor } from "./point_editor.mjs";
import { ModalEditor } from "./modal_editor.mjs";
import { DecayCurveEditor } from "./decay_curve_editor.mjs";
import { expandedSizeMeta } from "./size_meta.mjs";

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
  let digits = 3;
  if (descriptor.unit === "Hz") digits = value >= 1000 ? 0 : 1;
  else if (descriptor.unit === "dB" || descriptor.unit === "dB/oct") digits = 1;
  else if (descriptor.unit === "s") digits = value < .1 ? 3 : 2;
  else if (descriptor.unit === "x" || descriptor.unit === "oct") digits = 2;
  const suffix = descriptor.unit ? ` ${descriptor.unit}` : "";
  return `${Number(value).toFixed(digits)}${suffix}`;
}

const erb = frequency => 21.4 * Math.log10(1 + .00437 * frequency);
const inverseErb = rate => (10 ** (rate / 21.4) - 1) / .00437;

const ControlHelp = {
  model_level_db: "Pre-limiter synthesis level. Double-click to restore automatic reference matching.",
  direct_gain: "Amount of near-field contact heard directly before the cymbal body develops.",
  impact_tone_noise: "Balances pitched stick ping against broadband contact noise.",
  impact_width: "Scales contact duration. Short contacts are sharper; broad contacts suit softer implements.",
  bloom_level: "Audible gain from the dispersion return into the modal body. Zero removes bloom excitation completely.",
  bloom_nonlinearity: "Signal-dependent phase drive and delay excursion inside the dispersion loop.",
  bloom_development: "Feedback persistence of the bloom before its energy reaches the modal body.",
  bloom_diffusion: "Strength of the serial all-pass stages. Lower values sound nearer and more focused.",
  field_gain: "Overall level of the unified modal body, independent of relative anchor energies.",
  body_brightness: "Broad tilt applied to modal anchor energy around a 4 kHz pivot.",
  field_turbulence: "Global resolved-to-diffuse trajectory. Per-anchor widths scale this value.",
  field_packet_spread: "Maximum ERB-frequency spread of stochastic satellites around every anchor.",
  field_phase_bandwidth: "Rate of passive random phase decorrelation inside each modal packet.",
  field_exchange: "Amount of energy-preserving exchange between neighbouring modal states.",
};

function helpFor(key) {
  if (ControlHelp[key]) return ControlHelp[key];
  if (key.startsWith("body_decay_seconds_")) {
    return "T60 at this frequency knot; every unified modal packet shares the interpolated curve.";
  }
  if (key.startsWith("body_decay_frequency_")) {
    return "Frequency location of this shared T60 knot.";
  }
  if (key.startsWith("resolved_frequency_")) {
    return "Centre frequency of this constructive modal packet; this need not represent one measured physical eigenmode.";
  }
  if (key.startsWith("resolved_level_")) {
    return "Relative drive energy assigned to this packet. The exact minimum deactivates it.";
  }
  if (key.startsWith("resolved_turbulence_")) {
    return "Local multiplier for global turbulence: zero stays coherent, one follows the global setting, and two diffuses earlier.";
  }
  if (key.endsWith("radiation_enabled")) {
    return "Enable this static output-observation filter; it does not change stored body energy.";
  }
  if (key.includes("low_cut")) return "Static observation high-pass; it does not remove stored modal energy.";
  if (key.includes("high_cut")) return "Static observation low-pass; it does not shorten modal decay.";
  if (key.includes("colour_frequency")) return "Centre frequency of the observation colour bell.";
  if (key.includes("colour_gain")) return "Boost or cut of the observation colour bell.";
  if (key.includes("colour_q")) return "Bandwidth of the observation colour bell.";
  return "Double-click the control to restore its fitted default.";
}

export class FitControls {
  constructor({ descriptors, state, onChange, onLevelReset }) {
    this.descriptors = descriptors;
    this.byKey = new Map(descriptors.map(item => [item.key, item]));
    this.state = state;
    this.onChange = onChange;
    this.onLevelReset = onLevelReset;
  }

  descriptor(key) {
    const result = this.byKey.get(key);
    if (!result) throw new Error(`missing crash control ${key}`);
    return result;
  }

  value(key) {
    const descriptor = this.descriptor(key);
    return this.state.macros[descriptor.index];
  }

  setValue(key, value) {
    const descriptor = this.descriptor(key);
    this.state.macros[descriptor.index] = clamp(
      value, descriptor.minimum, descriptor.maximum,
    );
    this.onChange(key);
  }

  refresh(key) {
    const descriptor = this.descriptor(key);
    const row = document.querySelector(`[data-fit-key="${key}"]`);
    if (!row) return;
    row.querySelector("input").value = normalized(descriptor, this.value(key));
    row.querySelector("output").textContent = valueText(descriptor, this.value(key));
  }

  build() {
    document.querySelectorAll("[data-fit-controls]").forEach(
      element => element.replaceChildren(),
    );
    this.slider("model_level_db", "model-level", {
      reset: () => this.onLevelReset(),
    });
    this.sliders("impact-controls", [
      "direct_gain", "impact_tone_noise", "impact_width",
    ], {
      impact_tone_noise: ["ping", "noise"], impact_width: ["short", "broad"],
    });
    this.sliders("bloom-controls", [
      "bloom_level", "bloom_nonlinearity", "bloom_development",
      "bloom_diffusion",
    ], {
      bloom_level: ["off", "loud"],
      bloom_nonlinearity: ["linear", "nonlinear"],
      bloom_development: ["immediate", "slow"],
      bloom_diffusion: ["focused", "diffuse"],
    });
    this.sliders("body-controls", ["body_tone_wash"], {
      body_tone_wash: ["resolved", "wash"],
    });
    this.buildResolvedEditor();
    this.buildBodyModel();
    this.buildDenseControls();
    this.buildDenseWashEditor();
    this.buildDecayEditor();
    this.buildTurbulence();
    this.buildRadiation();
  }

  buildBodyModel() {
    this.sliders("field-turbulence-controls", [
      "field_gain", "body_brightness", "field_turbulence",
    ], {
      field_gain: ["quiet", "loud"],
      body_brightness: ["dark", "bright"],
      field_turbulence: ["resolved", "turbulent"],
    }, () => this.resolvedEditor?.refresh());
    this.sliders("field-advanced-controls", [
      "field_packet_spread", "field_phase_bandwidth", "field_exchange",
    ], {
      field_packet_spread: ["tight", "broad"],
      field_phase_bandwidth: ["coherent", "diffuse"],
      field_exchange: ["independent", "coupled"],
    }, () => this.resolvedEditor?.refresh());
  }

  setBodyMode(mode, notify = true) {
    const descriptor = this.descriptor("unified_body_enabled");
    this.state.macros[descriptor.index] = mode === "unified" ? 1 : 0;
    if (notify) this.onChange("body_ui_mode");
  }

  bodyMode() {
    return this.value("unified_body_enabled") >= .5 ? "unified" : "legacy";
  }

  applySizeMeta(position) {
    for (const { descriptor, value } of expandedSizeMeta(
      this.descriptors, position,
    )) {
      this.state.macros[descriptor.index] = value;
    }
    this.build();
    this.onChange("size_meta");
  }

  sliders(parentId, keys, labels = {}, afterInput) {
    for (const key of keys) this.slider(key, parentId, {
      labels: labels[key], afterInput,
    });
  }

  slider(key, parentId, options = {}) {
    const descriptor = this.descriptor(key);
    const row = document.createElement("label");
    row.className = "slider-row";
    row.dataset.fitKey = key;
    row.dataset.tooltip = helpFor(key);
    const title = document.createElement("span");
    title.textContent = descriptor.name;
    const input = document.createElement("input");
    input.type = "range"; input.min = 0; input.max = 1; input.step = 1 / 500;
    input.value = normalized(descriptor, this.value(key));
    const output = document.createElement("output");
    output.textContent = valueText(descriptor, this.value(key));
    input.oninput = () => {
      const raw = denormalized(descriptor, Number(input.value));
      this.setValue(key, options.coerce ? options.coerce(raw) : raw);
      output.textContent = valueText(descriptor, this.value(key));
      options.afterInput?.();
    };
    input.ondblclick = event => {
      event.preventDefault();
      if (options.reset) options.reset();
      else this.setValue(key, descriptor.defaultValue);
      input.value = normalized(descriptor, this.value(key));
      output.textContent = valueText(descriptor, this.value(key));
      options.afterInput?.();
    };
    row.append(title, input, output);
    if (options.labels) {
      const endpoints = document.createElement("span");
      endpoints.className = "slider-endpoints";
      endpoints.innerHTML = `<i>${options.labels[0]}</i><i>${options.labels[1]}</i>`;
      row.append(endpoints);
    }
    document.getElementById(parentId).append(row);
  }

  curveDescriptors(prefix) {
    const frequencies = this.descriptors.filter(item =>
      item.key.startsWith(`${prefix}_frequency_`));
    const levels = this.descriptors.filter(item =>
      item.key.startsWith(`${prefix}_level_`) ||
      item.key.startsWith(`${prefix}_seconds_`));
    return { frequencies, levels };
  }

  orderedFrequency(descriptors, index, value) {
    const lower = index === 0 ? descriptors[index].minimum :
      this.state.macros[descriptors[index - 1].index] * 1.03;
    const upper = index + 1 === descriptors.length ? descriptors[index].maximum :
      this.state.macros[descriptors[index + 1].index] / 1.03;
    return clamp(value, Math.min(lower, upper), Math.max(lower, upper));
  }

  setCurvePoint(curve, index, frequency, level, changeKey = "curve") {
    const frequencyDescriptor = curve.frequencies[index];
    const levelDescriptor = curve.levels[index];
    this.state.macros[frequencyDescriptor.index] = this.orderedFrequency(
      curve.frequencies, index, frequency,
    );
    this.state.macros[levelDescriptor.index] = clamp(
      level, levelDescriptor.minimum, levelDescriptor.maximum,
    );
    this.onChange(changeKey);
  }

  buildCurveInspector(parentId, titleText, curve, index, editor) {
    const parent = document.getElementById(parentId);
    parent.replaceChildren();
    const title = document.createElement("b");
    title.textContent = `${titleText} ${index + 1}`;
    parent.append(title);
    this.slider(curve.frequencies[index].key, parentId, {
      coerce: value => this.orderedFrequency(curve.frequencies, index, value),
      afterInput: () => editor.paint(),
    });
    this.slider(curve.levels[index].key, parentId, {
      afterInput: () => editor.paint(),
    });
  }

  buildResolvedEditor() {
    const curve = this.curveDescriptors("resolved");
    const turbulence = this.descriptors.filter(item =>
      item.key.startsWith("resolved_turbulence_"));
    const points = () => curve.frequencies.map((frequency, index) => {
      const level = this.state.macros[curve.levels[index].index];
      return {
        frequency: this.state.macros[frequency.index], level,
        turbulence: this.state.macros[turbulence[index].index],
        active: level > curve.levels[index].minimum + 1.e-3,
      };
    });
    const replace = (next, changeKey) => {
      next.forEach((point, index) => {
        this.state.macros[curve.frequencies[index].index] = clamp(
          point.frequency, curve.frequencies[index].minimum,
          curve.frequencies[index].maximum,
        );
        this.state.macros[curve.levels[index].index] = clamp(
          point.level, curve.levels[index].minimum, curve.levels[index].maximum,
        );
        this.state.macros[turbulence[index].index] = clamp(
          point.turbulence, turbulence[index].minimum, turbulence[index].maximum,
        );
      });
      this.onChange(changeKey);
    };
    let editor;
    const inspect = index => {
      const selection = document.getElementById("modal-selection");
      selection.replaceChildren();
      selection.classList.toggle("inactive", index === null);
      selection.setAttribute("aria-disabled", index === null ? "true" : "false");
      const inspected = index ?? 0;
      const title = document.createElement("b");
      title.textContent = index === null
        ? "No modal anchor selected"
        : `Selected modal anchor ${index + 1}`;
      selection.append(title);
      this.slider(curve.frequencies[inspected].key, "modal-selection", {
        coerce: value => editor.snapFrequency(value),
        afterInput: () => editor.refresh(),
      });
      this.slider(curve.levels[inspected].key, "modal-selection", {
        afterInput: () => editor.refresh(),
      });
      this.slider(turbulence[inspected].key, "modal-selection", {
        afterInput: () => editor.refresh(),
      });
      const names = ["Centre frequency", "Anchor level", "Local turbulence"];
      [...selection.querySelectorAll(".slider-row")].forEach((row, rowIndex) => {
        row.querySelector("span").textContent = names[rowIndex];
        if (index !== null) return;
        row.querySelector("input").disabled = true;
        row.querySelector("output").textContent = "—";
      });
    };
    editor = new ModalEditor(document.getElementById("modal-editor"), {
      minimumFrequency: 40, maximumFrequency: 15000,
      minimumLevel: -72, maximumLevel: 6, points,
      globalTurbulence: () => this.value("field_turbulence"),
      packetSpread: () => this.value("field_packet_spread"),
      replace,
      insert: (frequency, level) => {
        const next = points();
        const index = next.findIndex(point => !point.active);
        if (index < 0) return null;
        next[index] = {
          frequency, level: Math.max(level, -48), turbulence: 1, active: true,
        };
        replace(next, "modal_insert");
        return index;
      },
      remove: index => {
        const next = points();
        next[index].level = curve.levels[index].minimum;
        next[index].active = false;
        replace(next, "modal_delete");
      },
      select: inspect,
      readout: text => {
        document.getElementById("modal-readout").textContent = text;
      },
    });
    const tools = { edit: "edit", level: "shape", paint: "paint" };
    const setTool = tool => {
      editor.setTool(tool);
      for (const [name, value] of Object.entries(tools)) {
        const button = document.getElementById(`modal-tool-${name}`);
        button.classList.toggle("active", value === tool);
        button.setAttribute("aria-pressed", value === tool ? "true" : "false");
      }
    };
    document.getElementById("modal-tool-edit").onclick = () => setTool("edit");
    document.getElementById("modal-tool-level").onclick = () => setTool("shape");
    document.getElementById("modal-tool-paint").onclick = () => setTool("paint");
    const brush = document.getElementById("modal-brush");
    brush.oninput = () => {
      editor.setBrushWidth(brush.value);
      brush.nextElementSibling.textContent =
        `${Number(brush.value).toFixed(2)} ERB`;
    };
    document.getElementById("modal-clear").onclick = () => {
      const next = points();
      next.forEach((point, index) => {
        point.level = curve.levels[index].minimum;
        point.active = false;
      });
      replace(next, "modal_clear");
      editor.refresh(); inspect(null);
    };
    document.getElementById("modal-preset").onchange = event => {
      if (!event.target.value) return;
      this.applyModalPreset(event.target.value, curve, turbulence);
      event.target.value = "";
      editor.refresh(); inspect(null);
    };
    this.bindHarmonicGuide(editor);
    setTool("edit");
    this.resolvedCurve = curve;
    this.resolvedEditor = editor;
    inspect(null);
  }

  bindHarmonicGuide(editor) {
    const guide = document.getElementById("harmonic-guide");
    const note = document.getElementById("harmonic-note");
    const octave = document.getElementById("harmonic-octave");
    const snap = document.getElementById("harmonic-snap");
    const snapAll = document.getElementById("harmonic-snap-all");
    const output = document.getElementById("harmonic-frequency");
    const update = () => {
      const midi = 12 * (Number(octave.value) + 1) + Number(note.value);
      const fundamentalHz = 440 * 2 ** ((midi - 69) / 12);
      note.disabled = !guide.checked;
      octave.disabled = !guide.checked;
      snap.disabled = !guide.checked;
      snapAll.disabled = !guide.checked;
      output.textContent = `${fundamentalHz.toFixed(2)} Hz`;
      editor.setHarmonicGuide({
        visible: guide.checked, fundamentalHz,
        snap: guide.checked && snap.checked,
      });
    };
    guide.onchange = update;
    note.onchange = update;
    octave.onchange = update;
    snap.onchange = update;
    snapAll.onclick = () => editor.snapActiveToHarmonics();
    update();
  }

  applyModalPreset(name, curve, turbulence) {
    const count = curve.frequencies.length;
    const levels = Array(count).fill(-72);
    const frequencies = curve.frequencies.map(item => item.defaultValue);
    const widths = turbulence.map(item => item.defaultValue);
    if (name === "fitted") {
      curve.levels.forEach((item, index) => { levels[index] = item.defaultValue; });
    } else {
      const definitions = {
        even: {
          active: 18, low: 120, high: 14500, centre: .58,
          peak: -5, edge: -18,
        },
        low: {
          active: 14, low: 70, high: 8500, centre: .42,
          peak: -4, edge: -22,
        },
        shimmer: {
          active: 20, low: 350, high: 15000, centre: .76,
          peak: -3, edge: -20,
        },
      };
      const preset = definitions[name];
      if (!preset) return;
      const first = erb(preset.low); const last = erb(preset.high);
      for (let index = 0; index < preset.active; ++index) {
        const position = index / Math.max(1, preset.active - 1);
        frequencies[index] = inverseErb(first + position * (last - first));
        const distance = Math.min(
          1, Math.abs(position - preset.centre) /
            Math.max(preset.centre, 1 - preset.centre),
        );
        levels[index] = preset.peak + distance *
          (preset.edge - preset.peak);
        widths[index] = name === "low" ? .65 :
          name === "shimmer" ? 1.25 : 1;
      }
    }
    for (let index = 0; index < count; ++index) {
      this.state.macros[curve.frequencies[index].index] = frequencies[index];
      this.state.macros[curve.levels[index].index] = levels[index];
      this.state.macros[turbulence[index].index] = widths[index];
    }
    this.onChange(`modal_preset_${name}`);
  }

  buildDenseControls() {
    this.sliders("dense-controls", [
      "dense_minimum_frequency", "dense_maximum_frequency",
      "dense_mode_density",
    ], {
      dense_mode_density: ["sparse", "dense"],
    });
    this.sliders("dense-controls", [
      "dense_spacing_jitter", "dense_gain_spread", "dense_decay_spread",
    ], { dense_spacing_jitter: ["regular", "irregular"] });
  }

  buildDenseWashEditor() {
    const parent = document.getElementById("dense-wash-editor");
    const curve = this.curveDescriptors("dense_wash");
    const editor = new PointEditor(parent, {
      label: "Dense wash spectral colour",
      xMinimum: 40, xMaximum: 22000, xScale: "erb",
      yMinimum: -24, yMaximum: 24, connected: true,
      yTicks: [
        { value: -18, label: "-18" }, { value: 0, label: "0 dB" },
        { value: 18, label: "+18" },
      ],
      movableX: true,
      points: () => curve.frequencies.map((frequency, index) => ({
        x: this.state.macros[frequency.index],
        y: this.state.macros[curve.levels[index].index],
      })),
      setPoint: (index, frequency, level) => {
        this.setCurvePoint(curve, index, frequency, level, "dense_wash_colour");
        this.buildCurveInspector(
          "dense-wash-selection", "Selected colour knot", curve, index, editor,
        );
      },
      resetPoint: index => {
        this.state.macros[curve.frequencies[index].index] =
          curve.frequencies[index].defaultValue;
        this.state.macros[curve.levels[index].index] =
          curve.levels[index].defaultValue;
        this.onChange("dense_wash_colour");
        this.buildCurveInspector(
          "dense-wash-selection", "Selected colour knot", curve, index, editor,
        );
      },
      select: index => this.buildCurveInspector(
        "dense-wash-selection", "Selected colour knot", curve, index, editor,
      ),
    });
    this.denseWashCurve = curve;
    this.denseWashEditor = editor;
    this.buildCurveInspector(
      "dense-wash-selection", "Selected colour knot", curve, 0, editor,
    );
  }

  buildDecayEditor() {
    const parent = document.getElementById("decay-editor");
    const curve = this.curveDescriptors("body_decay");
    const active = this.descriptors.filter(item =>
      item.key.startsWith("body_decay_active_"));
    const nyquist = () => .5 * (this.state.reference?.sampleRate ?? 48000);
    const points = () => curve.frequencies.flatMap((frequency, slot) => {
      if (slot !== 0 && slot + 1 !== curve.frequencies.length &&
          this.state.macros[active[slot].index] < .5) return [];
      return [{
        slot,
        x: slot === 0 ? 0 : slot + 1 === curve.frequencies.length
          ? nyquist() : this.state.macros[frequency.index],
        y: Math.log2(this.state.macros[curve.levels[slot].index]),
        fixed: slot === 0 || slot + 1 === curve.frequencies.length,
      }];
    }).sort((left, right) => left.x - right.x);
    const replace = (next, changeKey = "body_decay") => {
      for (const point of next) {
        if (!point.fixed) {
          this.state.macros[curve.frequencies[point.slot].index] = point.x;
        }
        this.state.macros[curve.levels[point.slot].index] = clamp(
          2 ** point.y, curve.levels[point.slot].minimum,
          curve.levels[point.slot].maximum,
        );
      }
      this.onChange(changeKey);
    };
    let editor;
    const inspect = slot => this.buildDecayInspector(
      curve, active, slot, editor, nyquist,
    );
    editor = new DecayCurveEditor(parent, {
      nyquist,
      minimumLogSeconds: Math.log2(.02), maximumLogSeconds: Math.log2(20),
      yTicks: [
        { value: Math.log2(.1), label: ".1 s" },
        { value: 0, label: "1 s" },
        { value: Math.log2(10), label: "10 s" },
      ],
      points,
      setPoint: (slot, frequency, logSeconds) => {
        const next = points();
        const point = next.find(item => item.slot === slot);
        if (!point) return;
        point.x = point.fixed ? point.x : frequency;
        point.y = logSeconds;
        replace(next);
      },
      replace,
      insert: (frequency, logSeconds) => {
        if (points().length >= curve.frequencies.length) return null;
        const slot = active.findIndex((descriptor, index) =>
          index > 0 && index + 1 < active.length &&
          this.state.macros[descriptor.index] < .5);
        if (slot < 0) return null;
        this.state.macros[active[slot].index] = 1;
        this.state.macros[curve.frequencies[slot].index] = clamp(
          frequency, 40, Math.min(20000, nyquist() - 1),
        );
        this.state.macros[curve.levels[slot].index] = clamp(
          2 ** logSeconds, curve.levels[slot].minimum,
          curve.levels[slot].maximum,
        );
        this.onChange("body_decay_insert");
        return slot;
      },
      remove: slot => {
        if (slot === 0 || slot + 1 === active.length) return;
        this.state.macros[active[slot].index] = 0;
        this.onChange("body_decay_remove");
      },
      reset: slot => {
        this.state.macros[curve.levels[slot].index] =
          curve.levels[slot].defaultValue;
        this.onChange("body_decay_reset");
      },
      select: inspect,
      readout: text => {
        document.getElementById("decay-readout").textContent = text;
      },
    });
    this.decayCurve = curve;
    this.decayEditor = editor;
    inspect(0);
  }

  buildDecayInspector(curve, active, slot, editor, nyquist) {
    const parent = document.getElementById("decay-selection");
    parent.replaceChildren();
    const last = curve.frequencies.length - 1;
    const title = document.createElement("b");
    title.textContent = slot === 0 ? "DC boundary" : slot === last
      ? "Nyquist boundary" : `Selected T60 knot ${slot + 1}`;
    parent.append(title);
    this.slider(curve.frequencies[slot].key, "decay-selection", {
      coerce: value => editor.constrainFrequency(slot, value),
      afterInput: () => editor.paint(),
    });
    const frequencyRow = parent.querySelector(".slider-row");
    frequencyRow.querySelector("span").textContent = "Frequency";
    if (slot === 0 || slot === last) {
      frequencyRow.querySelector("input").disabled = true;
      frequencyRow.querySelector("output").textContent = slot === 0
        ? "DC" : `${(nyquist() / 1000).toFixed(2)} kHz`;
    }
    this.slider(curve.levels[slot].key, "decay-selection", {
      afterInput: () => editor.paint(),
    });
    parent.lastElementChild.querySelector("span").textContent = "T60";
    const remove = document.createElement("button");
    remove.type = "button";
    remove.disabled = slot === 0 || slot === last;
    remove.textContent = remove.disabled ? "Boundary knot" : "Delete knot";
    remove.title = remove.disabled
      ? "DC and Nyquist boundary knots cannot be deleted"
      : "Delete this interior knot (also: double-click it or press Delete)";
    remove.onclick = () => {
      this.state.macros[active[slot].index] = 0;
      this.onChange("body_decay_remove");
      editor.refresh();
    };
    parent.append(remove);
  }

  buildTurbulence() {
    const descriptor = this.descriptor("turbulence_enabled");
    const parent = document.getElementById("turbulence-toggle");
    const label = document.createElement("label");
    label.className = "checkbox-row";
    const input = document.createElement("input");
    input.type = "checkbox"; input.checked = this.value(descriptor.key) >= .5;
    const fields = document.getElementById("turbulence-fields");
    const update = () => {
      fields.disabled = !input.checked;
      this.setValue(descriptor.key, input.checked ? 1 : 0);
    };
    input.onchange = update;
    input.ondblclick = event => {
      event.preventDefault();
      input.checked = descriptor.defaultValue >= .5; update();
    };
    label.append(input, descriptor.name); parent.append(label);
    const curve = this.curveDescriptors("turbulence");
    const editor = new PointEditor(document.getElementById("turbulence-editor"), {
      label: "Turbulence spectral colour",
      xMinimum: 40, xMaximum: 20000,
      yMinimum: -18, yMaximum: 18, connected: true,
      yTicks: [
        { value: -12, label: "-12" }, { value: 0, label: "0 dB" },
        { value: 12, label: "+12" },
      ],
      movableX: true,
      points: () => curve.frequencies.map((frequency, index) => ({
        x: this.state.macros[frequency.index],
        y: this.state.macros[curve.levels[index].index],
      })),
      setPoint: (index, frequency, level) => {
        this.setCurvePoint(curve, index, frequency, level);
        this.buildCurveInspector(
          "turbulence-selection", "Selected colour knot", curve, index, editor,
        );
      },
      resetPoint: index => {
        this.state.macros[curve.frequencies[index].index] =
          curve.frequencies[index].defaultValue;
        this.state.macros[curve.levels[index].index] =
          curve.levels[index].defaultValue;
        this.onChange("turbulence_curve");
        this.buildCurveInspector(
          "turbulence-selection", "Selected colour knot", curve, index, editor,
        );
      },
      select: index => this.buildCurveInspector(
        "turbulence-selection", "Selected colour knot", curve, index, editor,
      ),
    });
    this.buildCurveInspector(
      "turbulence-selection", "Selected colour knot", curve, 0, editor,
    );
    this.sliders("turbulence-controls", [
      "turbulence_amount", "turbulence_persistence",
    ]);
    fields.disabled = !input.checked;
  }

  buildRadiation() {
    this.sliders("direct-radiation", ["direct_low_cut", "direct_high_cut"]);
    this.sliders("dense-radiation", [
      "dense_low_cut", "dense_high_cut",
      "dense_colour_frequency", "dense_colour_gain",
    ]);
    this.sliders("sparse-radiation", [
      "sparse_low_cut", "sparse_high_cut",
      "sparse_colour_frequency", "sparse_colour_gain",
    ]);
    for (const path of ["direct", "sparse", "dense"]) {
      this.checkbox(`${path}_radiation_enabled`, `${path}-radiation-advanced`);
      this.sliders(`${path}-radiation-advanced`, [
        `${path}_low_cut_q`, `${path}_colour_q`, `${path}_high_cut_q`,
      ]);
    }
    this.sliders("direct-radiation-advanced", [
      "direct_colour_frequency", "direct_colour_gain",
    ]);
  }

  checkbox(key, parentId, onToggle) {
    const descriptor = this.descriptor(key);
    const label = document.createElement("label");
    label.className = "checkbox-row";
    label.dataset.fitKey = key;
    const input = document.createElement("input");
    input.type = "checkbox"; input.checked = this.value(key) >= .5;
    input.onchange = () => {
      this.setValue(key, input.checked ? 1 : 0);
      onToggle?.(input.checked);
    };
    input.ondblclick = event => {
      event.preventDefault();
      input.checked = descriptor.defaultValue >= .5;
      this.setValue(key, descriptor.defaultValue);
      onToggle?.(input.checked);
    };
    label.append(input, descriptor.name);
    document.getElementById(parentId).append(label);
    onToggle?.(input.checked);
  }

}
