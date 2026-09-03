export class WasmCrashEngine {
  constructor(module, sampleRate, readDescriptors = true) {
    this.module = module;
    this.handle = 0;
    this.sampleRate = 0;
    this.scratchFrames = 512;
    this.scratch = module._malloc(
      this.scratchFrames * Float32Array.BYTES_PER_ELEMENT,
    );
    if (!this.scratch) throw new Error("could not allocate DSP scratch space");
    this.macros = readDescriptors ? this.#readMacroDescriptors() : [];
    this.setSampleRate(sampleRate);
  }

  #readMacroDescriptors() {
    const result = [];
    for (let index = 0; index < this.module._tf_crash_macro_count(); ++index) {
      result.push({
        index,
        key: this.module.UTF8ToString(this.module._tf_crash_macro_key(index)),
        name: this.module.UTF8ToString(this.module._tf_crash_macro_name(index)),
        unit: this.module.UTF8ToString(this.module._tf_crash_macro_unit(index)),
        scale: ["linear", "logarithmic", "boolean"][
          this.module._tf_crash_macro_scale(index)
        ],
        minimum: this.module._tf_crash_macro_minimum(index),
        maximum: this.module._tf_crash_macro_maximum(index),
        defaultValue: this.module._tf_crash_macro_default(index),
      });
    }
    return result;
  }

  setSampleRate(sampleRate) {
    if (this.sampleRate === sampleRate && this.handle) return;
    if (this.handle) this.module._tf_crash_destroy(this.handle);
    this.handle = this.module._tf_crash_create(sampleRate);
    if (!this.handle) throw new Error(`unsupported sample rate: ${sampleRate}`);
    this.sampleRate = sampleRate;
  }

  setMacros(values) {
    const count = this.module._tf_crash_macro_count();
    for (let index = 0; index < count; ++index) {
      const value = values[index] ?? this.module._tf_crash_macro_default(index);
      if (!this.module._tf_crash_macro_set(this.handle, index, value)) {
        throw new Error(`could not set crash control ${index}`);
      }
    }
    if (!this.module._tf_crash_macro_commit(this.handle)) {
      throw new Error("could not commit crash controls");
    }
  }

  trigger({
    strength = .8, location = 1, hardness = .65, implement = .75,
    contactSpread = .2, seed = 1,
  }) {
    if (!this.module._tf_crash_trigger(
      this.handle, strength, location, hardness, implement, contactSpread, seed,
    )) throw new Error("could not trigger crash DSP");
  }

  reset() {
    if (!this.module._tf_crash_reset(this.handle)) throw new Error("reset failed");
  }

  processTo(destination, offset, frames) {
    for (let first = 0; first < frames; first += this.scratchFrames) {
      const count = Math.min(this.scratchFrames, frames - first);
      this.#process(count);
      const start = this.scratch >>> 2;
      destination.set(
        this.module.HEAPF32.subarray(start, start + count), offset + first,
      );
    }
  }

  processQuantum(destination) {
    if (destination.length > this.scratchFrames) {
      throw new Error(`audio quantum exceeds ${this.scratchFrames} frames`);
    }
    this.#process(destination.length);
    const heap = this.module.HEAPF32;
    const start = this.scratch >>> 2;
    for (let frame = 0; frame < destination.length; ++frame) {
      destination[frame] = heap[start + frame];
    }
  }

  #process(frames) {
    if (!this.module._tf_crash_process(this.handle, this.scratch, frames)) {
      throw new Error("crash render failed");
    }
  }
}
