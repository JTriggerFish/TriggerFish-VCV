import createModule from "./triggerfish-percussion.mjs";
import { WasmCrashEngine } from "./wasm_engine_core.mjs";

export class CrashEngine extends WasmCrashEngine {
  static async create(sampleRate = 48000) {
    const module = await createModule();
    return new CrashEngine(module, sampleRate);
  }

  render({
    seconds, strength, location, hardness, implement, contactSpread, seed, macros,
  }) {
    this.setMacros(macros);
    this.trigger({
      strength, location, hardness, implement, contactSpread, seed,
    });

    const frames = Math.max(1, Math.round(seconds * this.sampleRate));
    const result = new Float32Array(frames);
    this.processTo(result, 0, frames);
    return result;
  }
}
