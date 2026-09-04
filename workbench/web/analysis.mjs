const TAU = 2 * Math.PI;

export const windows = Object.freeze(["hann", "blackman-harris", "rectangular"]);

export function isPowerOfTwo(value) {
  return value >= 2 && (value & (value - 1)) === 0;
}

export function windowSamples(name, size) {
  if (!isPowerOfTwo(size)) throw new Error("FFT size must be a power of two");
  const result = new Float64Array(size);
  const denominator = Math.max(1, size - 1);
  for (let index = 0; index < size; ++index) {
    const phase = TAU * index / denominator;
    if (name === "hann") result[index] = 0.5 - 0.5 * Math.cos(phase);
    else if (name === "blackman-harris") {
      result[index] = 0.35875 - 0.48829 * Math.cos(phase) +
        0.14128 * Math.cos(2 * phase) - 0.01168 * Math.cos(3 * phase);
    } else if (name === "rectangular") result[index] = 1;
    else throw new Error(`unknown window: ${name}`);
  }
  return result;
}

export function fft(real, imaginary) {
  const size = real.length;
  if (!isPowerOfTwo(size) || imaginary.length !== size) {
    throw new Error("FFT arrays need one equal power-of-two size");
  }
  for (let index = 1, reversed = 0; index < size; ++index) {
    let bit = size >> 1;
    while (reversed & bit) { reversed ^= bit; bit >>= 1; }
    reversed ^= bit;
    if (index < reversed) {
      [real[index], real[reversed]] = [real[reversed], real[index]];
      [imaginary[index], imaginary[reversed]] =
        [imaginary[reversed], imaginary[index]];
    }
  }
  for (let length = 2; length <= size; length <<= 1) {
    const angle = -TAU / length;
    const stepReal = Math.cos(angle);
    const stepImaginary = Math.sin(angle);
    for (let first = 0; first < size; first += length) {
      let rotationReal = 1;
      let rotationImaginary = 0;
      for (let offset = 0; offset < length / 2; ++offset) {
        const even = first + offset;
        const odd = even + length / 2;
        const oddReal = real[odd] * rotationReal -
          imaginary[odd] * rotationImaginary;
        const oddImaginary = real[odd] * rotationImaginary +
          imaginary[odd] * rotationReal;
        real[odd] = real[even] - oddReal;
        imaginary[odd] = imaginary[even] - oddImaginary;
        real[even] += oddReal;
        imaginary[even] += oddImaginary;
        const nextReal = rotationReal * stepReal -
          rotationImaginary * stepImaginary;
        rotationImaginary = rotationReal * stepImaginary +
          rotationImaginary * stepReal;
        rotationReal = nextReal;
      }
    }
  }
}

export function stft(samples, sampleRate, options = {}) {
  const size = options.size ?? 2048;
  const hop = options.hop ?? size / 4;
  const floorDb = options.floorDb ?? -140;
  if (!isPowerOfTwo(size) || hop < 1 || sampleRate < 1) {
    throw new Error("invalid STFT settings");
  }
  const window = windowSamples(options.window ?? "hann", size);
  const windowSum = window.reduce((sum, value) => sum + value, 0);
  const bins = size / 2 + 1;
  const frames = Math.max(1, 1 + Math.ceil((samples.length - 1) / hop));
  const values = new Float32Array(frames * bins);
  const real = new Float64Array(size);
  const imaginary = new Float64Array(size);
  let peakDb = floorDb;
  for (let frame = 0; frame < frames; ++frame) {
    const first = frame * hop - size / 2;
    real.fill(0);
    imaginary.fill(0);
    for (let index = 0; index < size; ++index) {
      const source = first + index;
      if (source >= 0 && source < samples.length) {
        real[index] = samples[source] * window[index];
      }
    }
    fft(real, imaginary);
    for (let bin = 0; bin < bins; ++bin) {
      const edge = bin === 0 || bin === size / 2;
      const scale = edge ? windowSum : windowSum / 2;
      const magnitude = Math.hypot(real[bin], imaginary[bin]) /
        Math.max(scale, 1e-30);
      const valueDb = Math.max(
        floorDb, 20 * Math.log10(Math.max(magnitude, 1e-30)),
      );
      values[frame * bins + bin] = valueDb;
      peakDb = Math.max(peakDb, valueDb);
    }
  }
  return { values, frames, bins, size, hop, sampleRate, floorDb, peakDb };
}

export function centeredErrorStatistics(values) {
  if (!values.length) return { mean: 0, rmse: 0 };
  const mean = values.reduce((sum, value) => sum + value, 0) / values.length;
  const rmse = Math.sqrt(values.reduce(
    (sum, value) => sum + (value - mean) ** 2, 0,
  ) / values.length);
  return { mean, rmse };
}
