import createModule from "./triggerfish-percussion.mjs";

export class CrashEngine {
  static async create(sampleRate = 48000) {
    const module = await createModule();
    return new CrashEngine(module, sampleRate);
  }

  constructor(module, sampleRate) {
    this.module = module;
    this.handle = 0;
    this.sampleRate = 0;
    this.scratchFrames = 512;
    this.scratch = this.module._malloc(
      this.scratchFrames * Float32Array.BYTES_PER_ELEMENT,
    );
    if (!this.scratch) throw new Error("could not allocate DSP scratch space");
    this.macros = this.#readMacroDescriptors();
    this.setSampleRate(sampleRate);
  }

  #readMacroDescriptors() {
    const result = [];
    for (let index = 0; index < this.module._tf_crash_macro_count(); ++index) {
      result.push({
        index,
        name: this.module.UTF8ToString(this.module._tf_crash_macro_name(index)),
        unit: this.module.UTF8ToString(this.module._tf_crash_macro_unit(index)),
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
    for (const descriptor of this.macros) {
      const value = values[descriptor.index] ?? descriptor.defaultValue;
      if (!this.module._tf_crash_macro_set(this.handle, descriptor.index, value)) {
        throw new Error(`could not set ${descriptor.name}`);
      }
    }
    if (!this.module._tf_crash_macro_commit(this.handle)) {
      throw new Error("could not commit crash controls");
    }
  }

  trigger({ strength, location, hardness, seed }) {
    if (!this.module._tf_crash_trigger(
      this.handle, strength, location, hardness, seed,
    )) throw new Error("could not trigger crash DSP");
  }

  reset() {
    if (!this.module._tf_crash_reset(this.handle)) throw new Error("reset failed");
  }

  processTo(destination, offset, frames) {
    for (let first = 0; first < frames; first += this.scratchFrames) {
      const count = Math.min(this.scratchFrames, frames - first);
      if (!this.module._tf_crash_process(this.handle, this.scratch, count)) {
        throw new Error("crash render failed");
      }
      const start = this.scratch >>> 2;
      destination.set(
        this.module.HEAPF32.subarray(start, start + count), offset + first,
      );
    }
  }

  render({ seconds, strength, location, hardness, seed, macros }) {
    this.setMacros(macros);
    this.trigger({ strength, location, hardness, seed });

    const frames = Math.max(1, Math.round(seconds * this.sampleRate));
    const result = new Float32Array(frames);
    this.processTo(result, 0, frames);
    return result;
  }
}
