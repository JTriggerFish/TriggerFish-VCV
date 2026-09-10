import {BloomTimingKeys} from "./bloom_timing_meta.mjs";
import {bloomRateNormalized} from "./bloom_control_scaling.mjs";

const Names = ["Diffusion strength", "Excitation tilt", "Excitation centre"];
const format = (key, value) => key === "body_excitation_centre"
  ? `${Math.round(value)} Hz` : key === "body_brightness"
    ? `${value.toFixed(1)} dB/oct` : value.toFixed(2);

export function createTimingView(parent, descriptors) {
  const button = document.createElement("button");
  button.type = "button"; button.className = "bloom-timing-button";
  button.textContent = "Bloom timing…";
  button.setAttribute("aria-expanded", "false");
  button.setAttribute("aria-haspopup", "dialog");
  const panel = document.createElement("div");
  panel.id = "bloom-timing-popup"; panel.className = "bloom-timing-popover";
  panel.popover = "auto"; panel.setAttribute("role", "dialog");
  panel.setAttribute("aria-labelledby", "bloom-timing-title");
  button.setAttribute("popovertarget", panel.id);
  button.setAttribute("aria-controls", panel.id);
  panel.innerHTML = `<header><strong id="bloom-timing-title">Bloom timing</strong>
    <button type="button" data-action="close" aria-label="Close bloom timing">×</button></header>
    <p>Move the bloom earlier or later by adjusting these three controls together.
      This is a relative guide, not an added delay.</p>
    <label class="bloom-timing-slider" data-fit-key="bloom_timing_meta">
      <span>Timing · meta</span><output>Centre</output>
      <input type="range" min="-1" max="1" step="0.002" value="0" aria-label="Bloom timing meta">
      <span class="slider-endpoints"><i>earlier</i><i>later</i></span>
    </label>
    <div class="bloom-timing-preview"></div>
    <div class="bloom-meta-status" role="status"></div>
    <div class="bloom-timing-actions"><button type="button" data-action="reset">Return to centre</button>
      <button type="button" data-action="centre">Set centre here</button></div>`;
  const rows = BloomTimingKeys.map((key, i) => {
    const row = document.createElement("div"); row.className = "bloom-preview-row";
    row.dataset.previewKey = key;
    row.innerHTML = `<span>${Names[i]}</span><output></output>
      <div class="bloom-preview-track" aria-hidden="true"><b></b><i></i></div>`;
    panel.querySelector(".bloom-timing-preview").append(row);
    return {key, row, descriptor:descriptors.find(d => d.key === key)};
  });
  parent.append(button, panel);
  const paint = (baseline, values) => rows.forEach(({key, row, descriptor:d}) => {
    row.querySelector("output").textContent = `${format(key, baseline[key])} → ${format(key, values[key])}`;
    const normalized = value => key === "bloom_rate" ? bloomRateNormalized(d, value) : d.scale === "logarithmic"
      ? Math.log(value / d.minimum) / Math.log(d.maximum / d.minimum)
      : (value - d.minimum) / (d.maximum - d.minimum);
    row.style.setProperty("--before", `${100 * normalized(baseline[key])}%`);
    row.style.setProperty("--after", `${100 * normalized(values[key])}%`);
  });
  return {button, panel, paint, input:panel.querySelector("input"),
    output:panel.querySelector(".bloom-timing-slider output"),
    status:panel.querySelector(".bloom-meta-status")};
}

export function placeTimingPopover(button, panel) {
  const anchor = button.getBoundingClientRect(), bounds = panel.getBoundingClientRect();
  const section = button.closest("#bloom-controls")?.getBoundingClientRect() ?? anchor;
  const margin = 12, gap = 10;
  const left = section.right + gap + bounds.width <= innerWidth - margin
    ? section.right + gap : Math.max(margin, Math.min(anchor.left, innerWidth - bounds.width - margin));
  const top = Math.max(margin, Math.min(anchor.top, innerHeight - bounds.height - margin));
  panel.style.left = `${left}px`; panel.style.top = `${top}px`;
}
