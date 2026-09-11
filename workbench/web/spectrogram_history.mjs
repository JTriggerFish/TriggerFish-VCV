// Display history only: never feed old pixels into analysis or saved audio.
export function overlaySpectrum(previous, next) {
  if (!next.incomplete || !previous?.frames) return next;
  const oldDuration = (previous.frames - 1) * previous.hop / previous.sampleRate;
  const sameGrid = previous.bins === next.bins && previous.size === next.size &&
    previous.hop === next.hop && previous.sampleRate === next.sampleRate;
  const frames = Math.max(next.frames, sameGrid ? previous.frames :
    1 + Math.floor(oldDuration * next.sampleRate / next.hop));
  const values = new Float32Array(frames * next.bins).fill(next.floorDb);
  if (sameGrid) values.set(previous.values.subarray(0, previous.frames * previous.bins));
  else {
    // Keep the old display correctly positioned if FFT resolution/rate changed.
    const bins = Array.from({length: next.bins}, (_, bin) =>
      Math.round(bin * next.sampleRate / next.size * previous.size / previous.sampleRate));
    for (let frame = next.frames; frame < frames; ++frame) {
      const old = Math.round(frame * next.hop / next.sampleRate * previous.sampleRate / previous.hop);
      if (old >= previous.frames) continue;
      for (let bin = 0; bin < next.bins; ++bin) {
        if (bins[bin] < previous.bins)
          values[frame * next.bins + bin] = previous.values[old * previous.bins + bins[bin]];
      }
    }
  }
  values.set(next.values.subarray(0, next.frames * next.bins));
  return {...next, values, frames, writeFrames: next.frames};
}

export function writeEdgeFraction(spectrum, viewport, offset = 0) {
  if (!spectrum?.incomplete) return null;
  const time = Math.max(0, (spectrum.writeFrames ?? spectrum.frames) - 1) *
    spectrum.hop / spectrum.sampleRate;
  const fraction = (time - viewport.start - offset) / (viewport.end - viewport.start);
  return fraction >= 0 && fraction <= 1 ? fraction : null;
}
