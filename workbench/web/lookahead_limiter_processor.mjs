import { limiterLookaheadSeconds } from "./limiter_config.mjs";
import { TruePeakDetector } from "./true_peak_detector.mjs";

const ceiling = 10 ** (-1 / 20);

class LookaheadLimiter extends AudioWorkletProcessor {
  constructor() {
    super();
    this.lookahead = Math.max(
      1, Math.round(limiterLookaheadSeconds * sampleRate));
    this.ringSize = this.lookahead + 1;
    this.release = Math.exp(-1 / (0.1 * sampleRate));
    this.rings = [];
    this.detectors = [];
    this.peaks = new Float64Array(this.ringSize + 1);
    this.indices = new Float64Array(this.ringSize + 1);
    this.head = 0;
    this.tail = 0;
    this.sample = 0;
    this.write = 0;
    this.gain = 1;
    this.gainDb = 0;
    this.hold = 0;
    this.meterCounter = 0;
    this.inputPeak = 0;
  }

  ensureChannels(count) {
    while (this.rings.length < count) {
      this.rings.push(new Float32Array(this.ringSize));
      this.detectors.push(new TruePeakDetector());
    }
  }

  updateGain(peak) {
    const count = this.peaks.length;
    while (this.head !== this.tail &&
        this.indices[this.head] < this.sample - this.lookahead)
      this.head = (this.head + 1) % count;
    while (this.head !== this.tail &&
        this.peaks[(this.tail + count - 1) % count] <= peak)
      this.tail = (this.tail + count - 1) % count;
    this.peaks[this.tail] = peak;
    this.indices[this.tail] = this.sample++;
    this.tail = (this.tail + 1) % count;
    const windowPeak = this.peaks[this.head];
    // Four-times peak sampling can miss a Nyquist sinusoid by cos(pi/8).
    // Reserve that headroom plus 0.1 dB FIR tolerance; no makeup gain. Hold
    // bridges the two peaks per 20 Hz cycle without adding audio latency.
    const protectedCeiling = ceiling * Math.cos(Math.PI / 8) * 10 ** (-0.1 / 20);
    const targetDb = Math.min(0, 20 * Math.log10(
      protectedCeiling / Math.max(windowPeak, 1e-30)));
    // Reconstructed alternating peaks differ by floating-point roundoff.
    // They must refresh hold too, not trigger recovery on every other cycle.
    if (targetDb <= this.gainDb + 1e-4) {
      this.gainDb = Math.min(this.gainDb, targetDb);
      this.hold = Math.ceil(0.03 * sampleRate);
    } else if (this.hold > 0) {
      --this.hold;
    } else {
      this.gainDb = Math.min(targetDb, this.release * this.gainDb);
    }
    this.gain = 10 ** (this.gainDb / 20);
  }

  process(inputs, outputs) {
    const input = inputs[0];
    const output = outputs[0];
    if (!output.length) return true;
    this.ensureChannels(output.length);
    for (let frame = 0; frame < output[0].length; ++frame) {
      let peak = 0;
      for (let channel = 0; channel < output.length; ++channel) {
        const source = input[Math.min(channel, input.length - 1)];
        const value = Number.isFinite(source?.[frame]) ? source[frame] : 0;
        this.rings[channel][this.write] = value;
        peak = Math.max(peak, this.detectors[channel].process(value));
      }
      this.updateGain(peak);
      this.inputPeak = Math.max(this.inputPeak, peak);
      const read = (this.write + 1) % this.ringSize;
      for (let channel = 0; channel < output.length; ++channel) {
        output[channel][frame] = this.rings[channel][read] * this.gain;
      }
      this.write = read;
      if (++this.meterCounter >= 2048) {
        this.port.postMessage({
          reductionDb: this.gainDb,
          inputPeakDb: 20 * Math.log10(Math.max(this.inputPeak, 1e-30)),
        });
        this.meterCounter = 0;
        this.inputPeak = 0;
      }
    }
    return true;
  }
}

registerProcessor("triggerfish-lookahead-limiter", LookaheadLimiter);
