export class WasmPercussionEngine {
  constructor(module, sampleRate, recipeIndex = 0, readDescriptors = true,
    preparedOnly = false) {
    this.module = module;
    this.handle = 0;
    this.sampleRate = 0;
    this.recipeIndex = recipeIndex;
    this.readDescriptors = readDescriptors;
    this.preparedOnly = preparedOnly;
    this.scratchFrames = 512;
    this.scratch = module._malloc(
      this.scratchFrames * Float32Array.BYTES_PER_ELEMENT,
    );
    if (!this.scratch) throw new Error("could not allocate DSP scratch space");
    try {
      this.#replaceHandle(recipeIndex, sampleRate);
    } catch (error) {
      module._free(this.scratch);
      this.scratch = 0;
      throw error;
    }
    this.parameters = readDescriptors ? this.#readParameterDescriptors() : [];
    this.macros = this.parameters;
    this.descriptorCache = new Map([[recipeIndex, this.parameters]]);
  }

  #readParameterDescriptors(handle = this.handle) {
    const result = [];
    const count = this.module._tf_percussion_parameter_count(handle);
    for (let index = 0; index < count; ++index) {
      result.push({
        index,
        key: this.module.UTF8ToString(
          this.module._tf_percussion_parameter_key(handle, index)),
        name: this.module.UTF8ToString(
          this.module._tf_percussion_parameter_name(handle, index)),
        unit: this.module.UTF8ToString(
          this.module._tf_percussion_parameter_unit(handle, index)),
        scale: ["linear", "logarithmic", "boolean", "choice"][
          this.module._tf_percussion_parameter_scale(handle, index)
        ],
        minimum: this.module._tf_percussion_parameter_minimum(
          handle, index),
        maximum: this.module._tf_percussion_parameter_maximum(
          handle, index),
        defaultValue: this.module._tf_percussion_parameter_default(
          handle, index),
      });
    }
    return result;
  }

  setSampleRate(sampleRate) {
    if (this.sampleRate === sampleRate && this.handle) return;
    this.#replaceHandle(this.recipeIndex, sampleRate);
  }

  setRecipe(recipeIndex) {
    if (recipeIndex === this.recipeIndex && this.handle) return;
    this.#replaceHandle(recipeIndex, this.sampleRate || 48000);
    if (this.readDescriptors) {
      this.parameters = this.#readParameterDescriptors();
      this.macros = this.parameters;
      this.descriptorCache.set(recipeIndex, this.parameters);
    }
  }

  descriptorsForRecipe(recipeIndex) {
    const cached = this.descriptorCache.get(recipeIndex);
    if (cached) return cached;
    const handle = this.module._tf_percussion_create(
      recipeIndex, this.sampleRate || 48000,
    );
    if (!handle) throw new Error(`recipe is unavailable: ${recipeIndex}`);
    try {
      const descriptors = this.#readParameterDescriptors(handle);
      this.descriptorCache.set(recipeIndex, descriptors);
      return descriptors;
    } finally {
      this.module._tf_percussion_destroy(handle);
    }
  }

  #replaceHandle(recipeIndex, sampleRate) {
    const create = this.preparedOnly
      ? this.module._tf_percussion_create_unprepared
      : this.module._tf_percussion_create;
    const replacement = create(recipeIndex, sampleRate);
    if (!replacement)
      throw new Error(`unsupported percussion recipe or sample rate: ${recipeIndex}`);
    if (this.handle) this.module._tf_percussion_destroy(this.handle);
    this.handle = replacement;
    this.recipeIndex = recipeIndex;
    this.sampleRate = sampleRate;
  }

  setParameters(values) {
    this.#requireConfigurable();
    this.#stageParameters(values);
    this.#commit();
  }

  setMacros(values) { this.setParameters(values); }

  setRouting(values) {
    this.#requireConfigurable();
    this.#stageRouting(values);
    this.#commit();
  }

  setConfiguration(values, routing) {
    this.#requireConfigurable();
    this.stageConfiguration(values, routing);
    this.#commit();
  }

  stageConfiguration(values, routing) {
    this.#requireConfigurable();
    this.#stageParameters(values);
    this.#stageRouting(routing);
  }

  #stageParameters(values) {
    const count = this.module._tf_percussion_parameter_count(this.handle);
    for (let index = 0; index < count; ++index) {
      const value = values[index] ??
        this.module._tf_percussion_parameter_default(this.handle, index);
      if (!this.module._tf_percussion_parameter_set(this.handle, index, value))
        throw new Error(`could not set percussion parameter ${index}`);
    }
  }

  #stageRouting(values = []) {
    const count = this.module._tf_percussion_route_count(this.handle);
    for (let index = 0; index < count; ++index) {
      if (!this.module._tf_percussion_route_set(
        this.handle, index, values[index] ?? 1,
      )) throw new Error(`could not set percussion route ${index}`);
    }
  }

  #commit() {
    if (!this.module._tf_percussion_commit(this.handle))
      throw new Error("could not commit percussion configuration");
  }

  #requireConfigurable() {
    if (this.preparedOnly)
      throw new Error("prepared-only renderer cannot compile parameters");
  }

  trigger({
    strength = .8, location = 1, hardness = .65, implement = .75,
    contactSpread = .2, constraint = 0, seed = 1,
  }) {
    if (!this.module._tf_percussion_set_mute(this.handle, constraint))
      throw new Error("could not set percussion constraint");
    if (!this.module._tf_percussion_trigger(
      this.handle, strength, location, hardness, implement, contactSpread, seed,
    )) throw new Error("could not trigger percussion DSP");
  }

  setConstraint(amount) {
    if (!this.module._tf_percussion_set_mute(this.handle, amount))
      throw new Error("could not set percussion constraint");
  }

  reset() {
    if (!this.module._tf_percussion_reset(this.handle))
      throw new Error("reset failed");
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
    if (destination.length > this.scratchFrames)
      throw new Error(`audio quantum exceeds ${this.scratchFrames} frames`);
    this.#process(destination.length);
    const heap = this.module.HEAPF32;
    const start = this.scratch >>> 2;
    for (let frame = 0; frame < destination.length; ++frame)
      destination[frame] = heap[start + frame];
  }

  #process(frames) {
    if (!this.module._tf_percussion_process(this.handle, this.scratch, frames))
      throw new Error("percussion render failed");
  }

  exportPrepared() {
    const byteSize = this.module._tf_percussion_prepared_size(this.handle);
    const pointer = this.module._malloc(byteSize);
    if (!byteSize || !pointer)
      throw new Error("could not allocate prepared DSP data");
    try {
      if (!this.module._tf_percussion_export_prepared(
        this.handle, pointer, byteSize,
      )) throw new Error("could not prepare percussion configuration");
      return this.module.HEAPU8.slice(pointer, pointer + byteSize);
    } finally {
      this.module._free(pointer);
    }
  }

  applyPrepared(bytes) {
    const source = bytes instanceof Uint8Array ? bytes : new Uint8Array(bytes);
    const pointer = this.module._malloc(source.byteLength);
    if (!source.byteLength || !pointer)
      throw new Error("could not allocate prepared DSP data");
    try {
      this.module.HEAPU8.set(source, pointer);
      if (!this.module._tf_percussion_apply_prepared(
        this.handle, pointer, source.byteLength,
      )) throw new Error("prepared percussion configuration was rejected");
    } finally {
      this.module._free(pointer);
    }
  }

  destroy() {
    if (this.handle) this.module._tf_percussion_destroy(this.handle);
    if (this.scratch) this.module._free(this.scratch);
    this.handle = 0;
    this.scratch = 0;
  }
}
