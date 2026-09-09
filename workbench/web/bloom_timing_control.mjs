import {BloomTimingKeys, bloomTimingValues} from "./bloom_timing_meta.mjs";

// The captured centre is UI gesture state, never an additional DSP parameter.
export function mountBloomTiming(parent, {read, apply, descriptors}) {
  const panel = document.createElement("div");
  panel.className = "bloom-timing-meta";
  panel.innerHTML = `<label class="slider-row" data-fit-key="bloom_timing_meta">
    <span>Bloom timing · meta</span>
    <input type="range" min="-1" max="1" step="0.002" value="0"
      aria-label="Bloom timing meta">
    <output>Centre</output>
    <span class="slider-endpoints"><i>earlier</i><i>later</i></span>
  </label><div class="bloom-meta-caption">
    <span>Moves excitation and diffusion. Hold decay can compensate damping.</span>
    <button type="button">Set centre</button>
  </div><div class="bloom-meta-status" role="status"></div>`;
  parent.append(panel);
  const input = panel.querySelector("input"), output = panel.querySelector("output");
  const status = panel.querySelector(".bloom-meta-status");
  let baseline;
  panel.querySelector("label").dataset.tooltip =
    "Relative to the captured patch. Later lowers the excitation shelf centre, darkens its slope and reduces diffusion strength. No delayed burst, gain compensation or exact millisecond guarantee. Double-click returns to the captured centre. Saved fits contain the resulting ordinary controls.";
  panel.querySelector("button").dataset.tooltip =
    "Use the current sound as the new middle of this timing gesture. Does not change the sound.";
  const rebase = () => {
    baseline = Object.fromEntries(BloomTimingKeys.map(key => [key, read(key)]));
    input.value = 0; output.textContent = "Centre";
    input.disabled = baseline.bloom_rate === 0;
    status.textContent = input.disabled ? "Enable diffusion strength below to use timing." :
      "Double-click to return · direct edits establish a new centre";
  };
  const change = () => {
    const position = Number(input.value);
    const {values, limited} = bloomTimingValues(baseline, position, descriptors);
    apply(values);
    output.textContent = position === 0 ? "Centre" :
      `${position > 0 ? "Later" : "Earlier"} ${Math.abs(position).toFixed(2)}`;
    status.textContent = limited.length ? `At control limit: ${limited.join(", ")}` :
      "Excitation and diffusion updated below · no level compensation";
  };
  input.oninput = change;
  input.ondblclick = event => { event.preventDefault(); input.value = 0; change(); };
  panel.querySelector("button").onclick = rebase;
  rebase();
  return {rebase};
}
