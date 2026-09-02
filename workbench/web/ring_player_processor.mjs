class RingPlayer extends AudioWorkletProcessor {
  constructor(options) {
    super();
    const buffers = options.processorOptions;
    this.audio = new Float32Array(buffers.audioBuffer);
    this.control = new Int32Array(buffers.controlBuffer);
    this.capacity = this.audio.length;
  }

  process(_inputs, outputs) {
    const output = outputs[0][0];
    let read = Atomics.load(this.control, 0) >>> 0;
    const write = Atomics.load(this.control, 1) >>> 0;
    const available = Math.min((write - read) >>> 0, output.length);
    for (let frame = 0; frame < available; ++frame) {
      output[frame] = this.audio[read++ % this.capacity];
    }
    output.fill(0, available);
    Atomics.store(this.control, 0, read | 0);
    Atomics.notify(this.control, 0, 1);
    if (available < output.length) Atomics.add(this.control, 2, 1);
    return true;
  }
}

registerProcessor("triggerfish-ring-player", RingPlayer);
