// Discrete character choice and a slow-rate-friendly display mapping.
export const hasPairedRing = (layout, rate) => [3, 4].includes(Math.round(layout)) && rate > 0;
export const beatRatePosition = (maximum, value) => Math.log1p(value / .1) / Math.log1p(maximum / .1);
export const beatRateValue = (maximum, position) => .1 * Math.expm1(position * Math.log1p(maximum / .1));
export const beatDepthPosition = value => Math.sqrt(Math.max(0, Math.min(1, value)));
export const beatDepthValue = position => Math.max(0, Math.min(1, position)) ** 2;
export const ringBeatRate = (frequency, rate, tilt) => Math.max(0,
  Math.min(80, rate * (Math.max(1, frequency) / 125) ** Math.max(-1, Math.min(1, tilt))));

export function mountPacketLayout(parent, {read, set, reset}) {
  const row = document.createElement("label");
  row.className = "slider-row";
  row.dataset.fitKey = "field_distribution";
  row.dataset.tooltip = "Paired ring adds a gentle pulse to each clear ring; set its speed with Beat rate. Scattered and Even coverage shape the surrounding shimmer. Beating doublets pairs that shimmer instead of the main ring. For a breathing low note and sizzling highs, try Paired ring with clearer lows and noisier highs.";
  row.innerHTML = '<span>Ring character</span><select aria-label="Ring character"><option value="3">Paired ring</option><option value="4">Plate cloud</option><option value="0">Scattered</option><option value="1">Even coverage</option><option value="2">Beating doublets</option></select>';
  row.dataset.tooltip += " Plate cloud spreads stable side modes approximately evenly in Hz. Try one wide upper handle with high Sideband allocation for a dense metallic layer; its centre sets the cloud's tuning.";
  const select = row.querySelector("select");
  select.value = Math.round(read());
  select.onchange = () => set(Number(select.value));
  select.ondblclick = () => { reset(); select.value = Math.round(read()); };
  parent.append(row);
  return () => { select.value = Math.round(read()); };
}
