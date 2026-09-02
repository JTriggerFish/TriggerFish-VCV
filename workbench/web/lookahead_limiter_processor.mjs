const ceiling = 10 ** (-1 / 20);

class LookaheadLimiter extends AudioWorkletProcessor {
  constructor() {
    super();
    this.lookahead = Math.max(1, Math.round(0.005 * sampleRate));
    this.ringSize = this.lookahead + 1;
    this.release = Math.exp(-1 / (0.1 * sampleRate));
    this.rings = [];
    this.write = 0;
    this.gain = 1;
    this.hold = 0;
    this.meterCounter = 0;
  }

  ensureChannels(count) {
    while (this.rings.length < count) {
      this.rings.push(new Float32Array(this.ringSize));
    }
  }

  updateGain(peak) {
    const target = peak > ceiling ? ceiling / peak : 1;
    if (target < this.gain) {
      this.gain = target;
      this.hold = this.lookahead;
    } else if (this.hold > 0) {
      --this.hold;
    } else {
      this.gain = this.release * this.gain + (1 - this.release);
    }
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
        peak = Math.max(peak, Math.abs(value));
      }
      this.updateGain(peak);
      const read = (this.write + 1) % this.ringSize;
      for (let channel = 0; channel < output.length; ++channel) {
        output[channel][frame] = this.rings[channel][read] * this.gain;
      }
      this.write = read;
      if (++this.meterCounter >= 2048) {
        this.port.postMessage(20 * Math.log10(Math.max(this.gain, 1e-9)));
        this.meterCounter = 0;
      }
    }
    return true;
  }
}

registerProcessor("triggerfish-lookahead-limiter", LookaheadLimiter);
