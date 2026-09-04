export function waveformEnvelope(
  samples, sampleRate, startSeconds = 0, bucketCount = 3000,
) {
  const step = Math.max(1, Math.ceil(samples.length / bucketCount));
  const x = [];
  const y = [];
  for (let index = 0; index < samples.length; index += step) {
    const end = Math.min(samples.length, index + step);
    let minimum = samples[index];
    let maximum = samples[index];
    let minimumIndex = index;
    let maximumIndex = index;
    for (let sample = index + 1; sample < end; ++sample) {
      if (samples[sample] < minimum) {
        minimum = samples[sample]; minimumIndex = sample;
      }
      if (samples[sample] > maximum) {
        maximum = samples[sample]; maximumIndex = sample;
      }
    }
    const append = sample => {
      x.push(sample / sampleRate - startSeconds);
      y.push(samples[sample]);
    };
    if (minimumIndex < maximumIndex) {
      append(minimumIndex); append(maximumIndex);
    } else {
      append(maximumIndex); append(minimumIndex);
    }
  }
  return { x, y };
}

export function drawWaveform(state) {
  const traces = [];
  if (state.reference) traces.push({
    ...waveformEnvelope(
      state.reference.samples, state.reference.sampleRate,
      state.reference.cell?.onset_seconds ?? 0,
    ),
    name: "Reference", mode: "lines", line: { color: "#e8b45a", width: 1 },
  });
  if (state.synthesis) traces.push({
    ...waveformEnvelope(
      state.synthesis, state.reference?.sampleRate ?? 48000,
    ),
    name: "Synthesis", mode: "lines", line: { color: "#68a7d8", width: 1 },
  });
  const referenceSeconds = state.reference
    ? state.reference.samples.length / state.reference.sampleRate -
      (state.reference.cell?.onset_seconds ?? 0)
    : 0;
  const synthesisSeconds = state.synthesis
    ? state.synthesis.length / (state.reference?.sampleRate ?? 48000) : 0;
  Plotly.react("waveform", traces, {
    margin: { l: 44, r: 12, t: 8, b: 28 }, paper_bgcolor: "#0d1015",
    plot_bgcolor: "#0d1015", font: { color: "#94a0af", size: 10 },
    xaxis: {
      title: "seconds", gridcolor: "#252d38", zeroline: false,
      fixedrange: true, range: [0, Math.max(referenceSeconds, synthesisSeconds)],
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
