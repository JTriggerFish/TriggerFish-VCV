import { ModalEditor } from "./modal_editor.mjs";
import { MiniEqEditor } from "./mini_eq_editor.mjs";
import { turbulenceIntensity } from "./turbulence_profile.mjs";
import { mountModalTemplates } from "./modal_template_controls.mjs";
import { DecayCurveEditor } from "./decay_curve_editor.mjs";
import { expandedSizeMeta } from "./size_meta.mjs";
import { mountBloomTiming } from "./bloom_timing_control.mjs";
import { BloomTimingKeys } from "./bloom_timing_meta.mjs";
import { bloomRateNormalized, bloomRateDenormalized } from "./bloom_control_scaling.mjs";
import { mountPacketLayout, hasPairedRing, beatRatePosition, beatRateValue, ringBeatRate, beatDepthPosition, beatDepthValue } from "./packet_layout_control.mjs";
import { packetAllocation } from "./packet_allocation.mjs";
import { helpFor } from "./fit_control_help.mjs";
import { mountDecayHold } from "./decay_hold_control.mjs";

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
  if (["bloom_energy_acceleration", "bloom_energy_sensitivity"].includes(descriptor.key)
      && value > 0 && value < .01) digits = 5;
  if (descriptor.unit === "Hz") digits = value >= 1000 ? 0 : value < 10 ? 2 : 1;
  else if (descriptor.unit === "dB" || descriptor.unit === "dB/oct") digits = 1;
  else if (descriptor.unit === "s") digits = value < .1 ? 3 : 2;
  else if (descriptor.unit === "x" || descriptor.unit === "oct") digits = 2;
  const suffix = descriptor.unit ? ` ${descriptor.unit}` : "";
  return `${Number(value).toFixed(digits)}${suffix}`;
}



export class FitControls {
  constructor({ descriptors, state, onChange, onLevelReset, decayHold }) {
    this.decayHold = decayHold;
    this.descriptors = descriptors;
    this.byKey = new Map(descriptors.map(item => [item.key, item]));
    this.state = state;
    this.onChange = onChange;
    this.onLevelReset = onLevelReset;
    this.refreshers = new Map();
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
    if (BloomTimingKeys.includes(key)) this.bloomTiming?.rebase();
    this.onChange(key);
    if (/^output_(low_cut|high_cut|colour_|eq_)/.test(key))
      this.refreshRadiation();
  }

  refresh(key) {
    if (this.refreshers.has(key)) return this.refreshers.get(key)();
    const descriptor = this.descriptor(key);
    const row = document.querySelector(`[data-fit-key="${key}"]`);
    if (!row) return;
    row.querySelector("input").value = normalized(descriptor, this.value(key));
    row.querySelector("output").textContent = valueText(descriptor, this.value(key));
  }

  build() {
    this.eqEditors = [];
    this.holdControl?.destroy();
    this.bloomTiming?.destroy();
    this.decayEditor?.destroy();
    this.resolvedEditor?.destroy();
    this.modalTemplates?.destroy();
    this.refreshers.clear();
    document.querySelectorAll("[data-fit-controls]").forEach(
      element => element.replaceChildren(),
    );
    this.slider("model_level_db", "model-level", {
      reset: () => this.onLevelReset(),
    });
    this.sliders("impact-controls", [
      "impact_tone_noise", "impact_width",
    ], {
      impact_tone_noise: ["ping", "noise"], impact_width: ["short", "broad"],
    });
    this.sliders("impact-advanced-controls", [
      "impact_chirp_pitch", "impact_noise_tilt", "impact_micro_density",
      "velocity_brightness",
    ]);
    this.buildBloom();
    this.buildResolvedEditor();
    this.buildBodyModel();
    this.buildDecayEditor();
    this.buildRadiation();
  }

  buildBloom() {
    if (this.decayHold) this.holdControl = mountDecayHold(
      document.getElementById("bloom-controls"), {
        state: this.state, ...this.decayHold,
        mount: document.getElementById("bloom-hold-controls"),
        apply: values => {
          const changed = this.descriptors.filter(d => values[d.index] !== this.state.macros[d.index]);
          this.state.macros.splice(0, this.state.macros.length, ...values);
          changed.forEach(d => this.refresh(d.key));
          this.decayEditor?.refresh();
          this.onChange("hold_decay");
        },
      });
    this.bloomTiming = mountBloomTiming(document.getElementById("bloom-timing-controls"), {
      read: key => this.value(key), descriptors: this.descriptors,
      apply: values => {
        for (const [key, value] of Object.entries(values)) {
          this.state.macros[this.descriptor(key).index] = value;
          this.refresh(key);
        }
        this.onChange("bloom_timing_meta");
      },
    });
    this.slider("bloom_rate", "bloom-diffusion-controls", {
      labels: ["off", "strong"],
      normalize: bloomRateNormalized, denormalize: bloomRateDenormalized,
    });
    this.slider("bloom_energy_acceleration", "bloom-diffusion-controls", {
      labels: ["even", "concentrated"],
      normalize: (descriptor, value) => Math.cbrt(value),
      denormalize: (descriptor, position) => position ** 3,
    });
    this.slider("bloom_energy_sensitivity", "bloom-diffusion-controls", {
      labels: ["independent", "strong"],
      normalize: (descriptor, value) => Math.sqrt(value / 2),
      denormalize: (descriptor, position) => 2 * position ** 2,
    });
    this.sliders("bloom-excitation-controls", ["body_brightness", "body_excitation_centre"], {
      body_brightness: ["dark", "bright"],
    });
  }

  buildBodyModel() {
    this.sliders("field-tuning-controls", ["body_excitation", "body_tune"]);
    this.buildPacketTexture();
    this.buildBeatingControls();
    this.buildPhaseMovement();
    this.updateDoubletControl();
  }

  buildPacketTexture() {
    this.slider("field_turbulence", "field-turbulence-controls", {
      normalize: (d,v) => Math.log1p(v/.01)/Math.log1p(d.maximum/.01),
      denormalize: (d,p) => .01*Math.expm1(p*Math.log1p(d.maximum/.01)),
      afterInput: () => this.resolvedEditor?.refresh(),
    });
    this.sliders("field-turbulence-controls", ["field_turbulence_slope"], {
      field_turbulence: ["tonal", "noisy"],
      field_turbulence_slope: ["noisy lows", "noisy highs"],
    }, () => this.resolvedEditor?.refresh());
    const refreshLayout = mountPacketLayout(document.getElementById("field-turbulence-controls"), {
      read: () => this.value("field_distribution"),
      set: value => { this.setValue("field_distribution", value); this.updateDoubletControl(); },
      reset: () => { this.setValue("field_distribution", this.descriptor("field_distribution").defaultValue); this.updateDoubletControl(); },
    });
    this.refreshers.set("field_distribution", refreshLayout);
    this.sliders("field-turbulence-controls", [
      "field_packet_spread", "field_satellite_density",
    ], {}, () => this.resolvedEditor?.refresh());
  }

  buildBeatingControls() {
    this.slider("field_doublet_split", "field-beating-controls", {
      labels: ["still / slow", "fast shimmer"],
      step: .0002,
      normalize: (d, v) => beatRatePosition(d.maximum, v),
      denormalize: (d, p) => beatRateValue(d.maximum, p),
      afterInput: () => this.resolvedEditor?.refresh(),
    });
    this.slider("field_beat_depth", "field-beating-controls", {
      labels: ["steady", "deep pulses"],
      normalize: (_descriptor, value) => beatDepthPosition(value),
      denormalize: (_descriptor, position) => beatDepthValue(position),
      afterInput: () => this.resolvedEditor?.refresh(),
    });
    this.sliders("field-beating-controls", ["field_beat_rate_tilt"], {},
      () => this.resolvedEditor?.refresh());
  }

  buildPhaseMovement() {
    this.slider("field_phase_bandwidth", "field-blur-controls", {
      labels: ["stable beating", "noise blur"],
      normalize: (descriptor, value) => Math.sqrt(value / descriptor.maximum),
      denormalize: (descriptor, position) => descriptor.maximum * position ** 2,
    });
    this.sliders("field-blur-controls", ["field_phase_tilt"], {
      field_phase_tilt: ["blur bass", "blur treble"],
    }, () => this.resolvedEditor?.refresh());
    this.slider("field_wander_hz", "field-drift-controls", {
      normalize: (d,v) => Math.log1p(v/.1)/Math.log1p(d.maximum/.1),
      denormalize: (d,p) => .1*Math.expm1(p*Math.log1p(d.maximum/.1)),
    });
    this.sliders("field-drift-controls", ["field_wander_rate"]);
  }

  destroy() {
    this.holdControl?.destroy();
    this.bloomTiming?.destroy();
    this.decayEditor?.destroy();
    this.resolvedEditor?.destroy?.();
    this.modalTemplates?.destroy();
  }

  updateDoubletControl() {
    const row = document.querySelector('[data-fit-key="field_doublet_split"]');
    if (!row) return;
    const layout = Math.round(this.value("field_distribution"));
    const inactive = layout !== 2 && layout !== 3;
    row.querySelector("input").disabled = inactive;
    row.style.opacity = inactive ? ".45" : "1";
    for (const key of ["field_beat_depth", "field_beat_rate_tilt"]) {
      const pairedRow = document.querySelector(`[data-fit-key="${key}"]`);
      if (!pairedRow) continue;
      pairedRow.querySelector("input").disabled = inactive;
      pairedRow.style.opacity = inactive ? ".45" : "1";
    }
    this.resolvedEditor?.refresh();
  }

  pairedRingReadout(point) {
    if (this.value("field_beat_depth") === 0) return " · steady ring";
    const frequency = point.frequency * this.value("body_tune");
    const rate = ringBeatRate(frequency, this.value("field_doublet_split"),
      this.value("field_beat_rate_tilt"));
    return ` · paired ring, ${rate.toFixed(2)} Hz`;
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
    row.dataset.tooltip = options.help ?? helpFor(key);
    const title = document.createElement("span");
    title.textContent = descriptor.name;
    const input = document.createElement("input");
    input.type = "range"; input.min = 0; input.max = 1; input.step = options.step ?? 1 / 500;
    const normalize = options.normalize ?? normalized;
    const denormalize = options.denormalize ?? denormalized;
    input.value = normalize(descriptor, this.value(key));
    const output = document.createElement("output");
    output.textContent = valueText(descriptor, this.value(key));
    this.refreshers.set(key, () => {
      input.value = normalize(descriptor, this.value(key));
      output.textContent = valueText(descriptor, this.value(key));
    });
    input.oninput = () => {
      const raw = denormalize(descriptor, Number(input.value));
      this.setValue(key, options.coerce ? options.coerce(raw) : raw);
      output.textContent = valueText(descriptor, this.value(key));
      options.afterInput?.();
    };
    input.ondblclick = event => {
      event.preventDefault();
      if (options.reset) options.reset();
      else this.setValue(key, descriptor.defaultValue);
      input.value = normalize(descriptor, this.value(key));
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

  buildResolvedEditor() {
    const curve = this.curveDescriptors("resolved");
    const minimumFrequency = curve.frequencies[0].minimum;
    const maximumFrequency = curve.frequencies[0].maximum;
    const turbulence = this.descriptors.filter(item =>
      item.key.startsWith("resolved_turbulence_"));
    const allocation = this.descriptors.filter(item =>
      item.key.startsWith("resolved_allocation_"));
    const points = () => curve.frequencies.map((frequency, index) => {
      const level = this.state.macros[curve.levels[index].index];
      return {
        frequency: this.state.macros[frequency.index], level,
        turbulence: this.state.macros[turbulence[index].index],
        allocation: this.state.macros[allocation[index].index],
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
        this.state.macros[allocation[index].index] = clamp(point.allocation ?? 1, 0, 4);
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
      this.slider(allocation[inspected].key, "modal-selection", {
        afterInput: () => editor.refresh(),
      });
      const names = ["Centre frequency", "Modal prominence", "Local noisiness", "Sideband allocation"];
      [...selection.querySelectorAll(".slider-row")].forEach((row, rowIndex) => {
        row.querySelector("span").textContent = names[rowIndex];
        if (index !== null) return;
        row.querySelector("input").disabled = true;
        row.querySelector("output").textContent = "—";
      });
    };
    editor = new ModalEditor(document.getElementById("modal-editor"), {
      minimumFrequency: Math.max(20, minimumFrequency), maximumFrequency, frequencyScale: "log",
      minimumLevel: -72, maximumLevel: 6, points,
      globalTurbulence: () => this.value("field_turbulence"),
      spectralTurbulence: frequency => turbulenceIntensity(frequency,
        this.value("field_turbulence"), this.value("field_turbulence_slope"),
        1000, 1,
        true),
      effectiveTurbulence: point => turbulenceIntensity(point.frequency,
        this.value("field_turbulence"), this.value("field_turbulence_slope"),
        1000, point.turbulence,
        true),
      packetSpread: () => this.value("field_packet_spread"),
      replace,
      insert: (frequency, level) => {
        const next = points();
        const index = next.findIndex(point => !point.active);
        if (index < 0) return null;
        next[index] = {
          frequency, level: Math.max(level, curve.levels[index].minimum + .1),
          turbulence: 1, active: true,
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
        const pool = packetAllocation(points(), this.value("field_satellite_density"),
          point => turbulenceIntensity(point.frequency, this.value("field_turbulence"),
            this.value("field_turbulence_slope"), 1000,
            point.turbulence, true) * this.value("field_packet_spread"),
          hasPairedRing(this.value("field_distribution"), this.value("field_doublet_split")));
        document.getElementById("modal-readout").textContent = `${text} · ${pool.count}/512 oscillators`;
        const title = document.querySelector("#modal-selection > b");
        if (title && editor && editor.selected !== null) {
          const count = pool.centresPerHandle + 2 * pool.pairs[editor.selected];
          const ring = pool.centresPerHandle === 2
            ? this.pairedRingReadout(points()[editor.selected]) : "";
          title.textContent = `Selected modal anchor ${editor.selected + 1} · ≈${count} ${count === 1 ? "oscillator" : "oscillators"}${ring}`;
        }
      },
    });
    // Preserve the shared guide controls across generator/control rebuilds.
    const guideToolbar = document.querySelector(".harmonic-toolbar");
    guideToolbar.remove();
    const templates = this.modalTemplates = mountModalTemplates(document.getElementById("modal-templates"), {
      capacity: curve.frequencies.length, minimumFrequency, maximumFrequency,
      defaultFamily: "harmonic", open: true, noisiness: true,
      apply: generated => {
        const next = points().map((point, i) => generated[i]
          ? {...point, ...generated[i], allocation:1} : {...point, level:-72, active:false});
        replace(next, "modal_template"); editor.select(null); editor.refresh();
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
      editor.select(null); editor.refresh();
    };
    document.querySelector("#modal-templates .template-panel").append(guideToolbar);
    this.bindHarmonicGuide(editor, templates, minimumFrequency);
    setTool("edit");
    this.resolvedCurve = curve;
    this.resolvedEditor = editor;
    inspect(null);
  }

  bindHarmonicGuide(editor, templates, minimumFrequency) {
    const guide = document.getElementById("harmonic-guide");
    const snap = document.getElementById("harmonic-snap");
    const snapAll = document.getElementById("harmonic-snap-all");
    const output = document.getElementById("harmonic-frequency");
    const update = () => {
      const fundamentalHz = templates.fundamentalHz;
      if (!Number.isFinite(fundamentalHz) || fundamentalHz < minimumFrequency ||
          fundamentalHz > editor.options.maximumFrequency) return;
      snap.disabled = !guide.checked;
      snapAll.disabled = !guide.checked;
      output.textContent = `${fundamentalHz.toFixed(2)} Hz`;
      editor.setHarmonicGuide({
        visible: guide.checked, fundamentalHz,
        snap: guide.checked && snap.checked,
      });
    };
    guide.onchange = update;
    templates.onPitchChange(update);
    snap.onchange = update;
    snapAll.onclick = () => editor.snapActiveToHarmonics();
    update();
  }


  buildDecayEditor() {
    const parent = document.getElementById("decay-editor");
    const interiorFrequencies = this.descriptors.filter(item =>
      item.key.startsWith("body_decay_frequency_"));
    const levels = this.descriptors.filter(item =>
      item.key.startsWith("body_decay_seconds_"));
    const interiorActive = this.descriptors.filter(item =>
      item.key.startsWith("body_decay_active_"));
    const curve = {
      frequencies: [null, ...interiorFrequencies, null], levels,
    };
    const active = [null, ...interiorActive, null];
    const decayMinimum = 40;
    const decayMaximum = 15000;
    const points = () => curve.frequencies.flatMap((frequency, slot) => {
      if (slot !== 0 && slot + 1 !== curve.frequencies.length &&
          this.state.macros[active[slot].index] < .5) return [];
      return [{
        slot,
        x: slot === 0 ? decayMinimum :
          slot + 1 === curve.frequencies.length
            ? decayMaximum : this.state.macros[frequency.index],
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
      curve, active, slot, editor, decayMinimum, decayMaximum,
    );
    editor = new DecayCurveEditor(parent, {
      minimumFrequency: decayMinimum,
      maximumFrequency: decayMaximum,
      minimumLogSeconds: Math.log2(levels[0].minimum),
      maximumLogSeconds: Math.log2(levels[0].maximum),
      yTicks: [
        { value: Math.log2(.1), label: ".1 s" },
        { value: 0, label: "1 s" },
        { value: Math.log2(10), label: "10 s" },
        { value: Math.log2(levels[0].maximum), label: `${levels[0].maximum} s` },
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
        const interior = interiorActive.findIndex(descriptor =>
          this.state.macros[descriptor.index] < .5);
        if (interior < 0) return null;
        const slot = interior + 1;
        this.state.macros[active[slot].index] = 1;
        this.state.macros[curve.frequencies[slot].index] = clamp(
          frequency, decayMinimum, decayMaximum,
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

  buildDecayInspector(
    curve, active, slot, editor, decayMinimum, decayMaximum,
  ) {
    const parent = document.getElementById("decay-selection");
    parent.replaceChildren();
    const last = curve.frequencies.length - 1;
    const title = document.createElement("b");
    title.textContent = slot === 0 ? "Low boundary" : slot === last
      ? "High boundary" : `Selected T60 knot ${slot + 1}`;
    parent.append(title);
    if (slot === 0 || slot === last) {
      const frequencyRow = document.createElement("div");
      frequencyRow.className = "slider-row fixed-readout";
      const label = document.createElement("span");
      label.textContent = "Frequency";
      const output = document.createElement("output");
      const boundary = slot === 0 ? decayMinimum : decayMaximum;
      output.textContent = boundary >= 1000
        ? `${(boundary / 1000).toFixed(1)} kHz` : `${boundary} Hz`;
      frequencyRow.append(label, output);
      parent.append(frequencyRow);
    } else {
      this.slider(curve.frequencies[slot].key, "decay-selection", {
        coerce: value => editor.constrainFrequency(slot, value),
        afterInput: () => editor.paint(),
      });
      parent.lastElementChild.querySelector("span").textContent = "Frequency";
    }
    this.slider(curve.levels[slot].key, "decay-selection", {
      afterInput: () => editor.paint(),
    });
    parent.lastElementChild.querySelector("span").textContent = "T60";
    const remove = document.createElement("button");
    remove.type = "button";
    remove.disabled = slot === 0 || slot === last;
    remove.textContent = "Delete knot";
    remove.title = remove.disabled
      ? "Low and high boundary knots cannot be deleted"
      : "Delete this interior knot (also: double-click it or press Delete)";
    remove.onclick = () => {
      this.state.macros[active[slot].index] = 0;
      this.onChange("body_decay_remove");
      editor.refresh();
    };
    parent.append(remove);
  }

  buildRadiation() {
    this.sliders("observation-mix", ["direct_gain", "field_gain"]);
    this.checkbox("output_eq_enabled", "output-eq-controls");
    const editor = new MiniEqEditor(document.getElementById("output-eq-editor"), {
      prefix:"output", state:this.state, read:key=>this.value(key),
      set:(key,value)=>this.setValue(key,value), descriptor:key=>this.descriptor(key),
    });
    this.eqEditors.push(editor);
    for(const key of editor.keys)this.refreshers.set(key,()=>editor.refresh());
  }

  refreshRadiation() {
    this.eqEditors?.forEach(editor=>editor.refresh());
  }

  checkbox(key, parentId, onToggle) {
    const descriptor = this.descriptor(key);
    const label = document.createElement("label");
    label.className = "checkbox-row";
    label.dataset.fitKey = key;
    label.dataset.tooltip = helpFor(key);
    const input = document.createElement("input");
    input.type = "checkbox"; input.checked = this.value(key) >= .5;
    this.refreshers.set(key, () => { input.checked = this.value(key) >= .5; });
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
    onToggle?.(input.checked, true);
  }

}
