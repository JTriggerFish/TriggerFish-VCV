import {modalTemplate, modalTemplateLimit} from "./modal_templates.mjs";
import {mountTemplateNote} from "./modal_template_note.mjs";
import {templateField} from "./modal_template_fields.mjs";
import {reportError} from "./error_banner.mjs";

export function mountModalTemplates(parent, options) {
  const {capacity, minimumFrequency, maximumFrequency, apply,
    defaultFamily = "membrane", open = false, noisiness = false} = options;
  parent.innerHTML = `<details class="template-panel"><summary>Generate modes</summary>
    <div class="template-pitch"></div><div class="template-shape"></div>
    <div class="template-actions"><button type="button">Replace modes</button>
    <span class="template-status" role="status" aria-live="polite"></span></div>
    <p class="control-help">Preview settings here, then replace. Existing damping and bloom stay unchanged.</p></details>`;
  const details = parent.firstElementChild; details.open = open;
  const pitch = parent.querySelector(".template-pitch");
  const shape = parent.querySelector(".template-shape");
  const formula = document.createElement("label"); formula.textContent = "Series ";
  const family = document.createElement("select");
  family.setAttribute("aria-label", "Modal formula");
  family.append(new Option("Harmonic", "harmonic"), new Option("Membrane", "membrane"));
  family.value = defaultFamily; formula.append(family); pitch.append(formula);
  const fields = {
    fundamental: templateField(pitch, "fundamental", "Hz", 55, minimumFrequency, maximumFrequency, "any", false),
    count: templateField(shape, "count", "Mode count", Math.min(16, capacity), 1, capacity, 1),
    stretch: templateField(shape, "stretch", "Upper-mode stretch", 0, 0, 1, .01),
    harmonicCore: templateField(shape, "harmonicCore", "Harmonic core · modes", 4, 1, 8, 1),
    rolloff: templateField(shape, "rolloff", "Falloff · dB/oct", 6, -12, 24, .5),
    level: templateField(shape, "level", "Top level · dB", 0, -60, 6, .5),
  };
  fields.stretch.parentElement.dataset.tooltip = "0 keeps the original series. The low harmonic core stays unchanged; higher modes bend progressively upward. Adding more modes does not retune existing ones. Writes editable frequencies, not a runtime rule.";
  fields.harmonicCore.parentElement.dataset.tooltip = "Number of low modes kept at the original series frequencies (four by default). Stretch begins smoothly above this core. For Membrane, these retain the membrane ratios rather than becoming harmonics.";
  if (noisiness) {
    fields.turbulence = templateField(shape, "turbulence", "Noisiness response", 1, 0, 2, .05);
    fields.turbulence.parentElement.dataset.tooltip = "Writes each mode's local noisiness multiplier. 1 follows the body noisiness controls; 0 makes centres pure and disables that response. Not frequency detuning.";
  }
  mountTemplateNote(pitch, fields.fundamental);
  const button = parent.querySelector("button");
  const status = parent.querySelector(".template-status");
  const preview = () => {
    const points = previewTemplate({fields, family, options, status, button});
    parent.dispatchEvent(new Event("template-pitch-change"));
    return points;
  };
  pitch.addEventListener("input", preview); shape.addEventListener("input", preview);
  family.onchange = preview;
  button.onclick = () => {
    const points = preview(); if (!points) return;
    try { apply(points); status.textContent = `Applied ${points.length} modes. All remain editable.`; }
    catch (error) {
      status.dataset.error = "true"; status.textContent = error.message;
      reportError(error, "Mode generator");
    }
  };
  // Report committed invalid entries, not every intermediate keystroke.
  parent.addEventListener("change", () => {
    preview();
    if (status.dataset.error === "true") reportError(status.textContent, "Mode generator");
  });
  preview();
  return {get fundamentalHz() { return Number(fields.fundamental.value); },
    onPitchChange: callback => parent.addEventListener("template-pitch-change", callback),
    open: name => {
    family.value = name; details.open = true; preview();
    details.scrollIntoView({block:"nearest"});
  }};
}

function previewTemplate({fields, family, options, status, button}) {
  const values = Object.fromEntries(Object.entries(fields).map(([key, input]) => [key, Number(input.value)]));
  const limit = modalTemplateLimit({...options, family:family.value,
    fundamental:values.fundamental, stretch:values.stretch, harmonicCore:values.harmonicCore});
  fields.count.setMaximum(limit);
  for (const input of Object.values(fields)) input.setAttribute("aria-invalid", String(!input.checkValidity()));
  status.dataset.error = "false";
  try {
    if (!limit) throw Error(`Check base pitch (${options.minimumFrequency}–${options.maximumFrequency} Hz), stretch (0–1) and harmonic core (1–8).`);
    if (!Number.isInteger(values.count) || values.count < 1 || values.count > limit)
      throw Error(`Choose 1–${limit} modes at this pitch/stretch (${family.value === "membrane" ? "16-root membrane formula" : "harmonic series"}; ${options.capacity} handles maximum).`);
    for (const input of Object.values(fields)) {
      const invalid = !input.checkValidity(); input.setAttribute("aria-invalid", String(invalid));
      if (invalid) throw Error(`${input.getAttribute("aria-label")}: enter ${input.min}–${input.max}.`);
    }
    const points = modalTemplate({...options, ...values, family:family.value});
    const off = points.filter(p => !p.active).length;
    const limited = points.filter(p => p.level === 6).length;
    status.textContent = `${points.length} modes · ${Math.round(points[0].frequency)}–${Math.round(points.at(-1).frequency)} Hz · max ${limit}` +
      (off ? ` · ${off} below level floor (off)` : "") + (limited ? ` · ${limited} at +6 dB ceiling` : "");
    button.disabled = false; return points;
  } catch (error) {
    button.disabled = true; status.dataset.error = "true"; status.textContent = error.message; return null;
  }
}
