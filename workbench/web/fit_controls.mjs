import { PointEditor } from "./point_editor.mjs";
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
    this.sliders("bloom-controls", ["bloom_amount", "bloom_development"], {
      bloom_amount: ["subtle", "strong"],
      bloom_development: ["immediate", "slow"],
    });
    this.sliders("body-controls", ["body_tone_wash", "body_brightness"], {
      body_tone_wash: ["resolved", "wash"],
      body_brightness: ["dark", "bright"],
    });
    this.buildResolvedEditor();
    this.buildDenseControls();
    this.buildDenseWashEditor();
    this.buildDecayEditor();
    this.buildTurbulence();
    this.buildRadiation();
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
    const parent = document.getElementById("resolved-editor");
    const curve = this.curveDescriptors("resolved");
    const editor = new PointEditor(parent, {
      label: "Resolved modal ridges",
      xMinimum: 40, xMaximum: 22000, xScale: "erb",
      yMinimum: -24, yMaximum: 24, bars: true, paint: true,
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
        this.setCurvePoint(curve, index, frequency, level, "resolved_modes");
        this.buildCurveInspector(
          "resolved-selection", "Selected mode", curve, index, editor,
        );
      },
      resetPoint: index => {
        this.state.macros[curve.frequencies[index].index] =
          curve.frequencies[index].defaultValue;
        this.state.macros[curve.levels[index].index] =
          curve.levels[index].defaultValue;
        this.onChange("resolved_modes");
        this.buildCurveInspector(
          "resolved-selection", "Selected mode", curve, index, editor,
        );
      },
      select: index => this.buildCurveInspector(
        "resolved-selection", "Selected mode", curve, index, editor,
      ),
    });
    this.resolvedCurve = curve;
    this.resolvedEditor = editor;
    this.buildCurveInspector(
      "resolved-selection", "Selected mode", curve, 0, editor,
    );
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
    const editor = new PointEditor(parent, {
      label: "Shared frequency-dependent body T60",
      xMinimum: 40, xMaximum: 20000,
      yMinimum: Math.log2(.02), yMaximum: Math.log2(20), connected: true,
      yTicks: [
        { value: Math.log2(.1), label: ".1 s" },
        { value: 0, label: "1 s" },
        { value: Math.log2(10), label: "10 s" },
      ],
      movableX: true,
      points: () => curve.frequencies.map((frequency, index) => ({
        x: this.state.macros[frequency.index],
        y: Math.log2(this.state.macros[curve.levels[index].index]),
      })),
      setPoint: (index, frequency, logSeconds) => {
        this.setCurvePoint(curve, index, frequency, 2 ** logSeconds);
        this.buildCurveInspector(
          "decay-selection", "Selected T60 knot", curve, index, editor,
        );
      },
      resetPoint: index => {
        this.state.macros[curve.frequencies[index].index] =
          curve.frequencies[index].defaultValue;
        this.state.macros[curve.levels[index].index] =
          curve.levels[index].defaultValue;
        this.onChange("body_decay");
        this.buildCurveInspector(
          "decay-selection", "Selected T60 knot", curve, index, editor,
        );
      },
      select: index => this.buildCurveInspector(
        "decay-selection", "Selected T60 knot", curve, index, editor,
      ),
    });
    this.decayCurve = curve;
    this.decayEditor = editor;
    this.buildCurveInspector(
      "decay-selection", "Selected T60 knot", curve, 0, editor,
    );
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
    for (const path of ["direct", "sparse", "dense"]) {
      this.checkbox(`${path}_radiation_enabled`, `${path}-radiation`);
      this.sliders(`${path}-radiation`, [
        `${path}_low_cut`, `${path}_low_cut_q`,
        `${path}_colour_frequency`, `${path}_colour_gain`, `${path}_colour_q`,
        `${path}_high_cut`, `${path}_high_cut_q`,
      ]);
    }
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
