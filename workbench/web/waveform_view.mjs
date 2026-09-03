function decimate(samples, sampleRate, limit = 5000) {
  const step = Math.max(1, Math.floor(samples.length / limit));
  const x = [];
  const y = [];
  for (let index = 0; index < samples.length; index += step) {
    x.push(index / sampleRate);
    y.push(samples[index]);
  }
  return { x, y };
}

export function drawWaveform(state) {
  const traces = [];
  if (state.reference) traces.push({
    ...decimate(state.reference.samples, state.reference.sampleRate),
    name: "Reference", mode: "lines", line: { color: "#e8b45a", width: 1 },
  });
  if (state.synthesis) traces.push({
    ...decimate(state.synthesis, state.reference?.sampleRate ?? 48000),
    name: "Synthesis", mode: "lines", line: { color: "#68a7d8", width: 1 },
  });
  Plotly.react("waveform", traces, {
    margin: { l: 44, r: 12, t: 8, b: 28 }, paper_bgcolor: "#0d1015",
    plot_bgcolor: "#0d1015", font: { color: "#94a0af", size: 10 },
    xaxis: {
      title: "seconds", gridcolor: "#252d38", zeroline: false,
      fixedrange: true,
    },
    yaxis: {
      title: "amplitude", gridcolor: "#252d38",
      zerolinecolor: "#46515e", fixedrange: true,
    },
    dragmode: false,
    legend: { orientation: "h", x: .75, y: 1 }, uirevision: "waveform-v1",
  }, {
    responsive: true, displaylogo: false, displayModeBar: false,
    scrollZoom: false,
  });
}
