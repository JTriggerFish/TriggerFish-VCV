import createModule from "./triggerfish-percussion-audio.mjs";
import { WasmCrashEngine } from "./wasm_engine_core.mjs";

// MODULARIZE consistently exposes a promise even when Wasm compilation itself is
// synchronous. Register the processor immediately, then attach each node to the
// shared module as soon as that promise settles.
const modulePromise = Promise.resolve(createModule());

function nowMs() {
  return globalThis.performance?.now() ?? currentTime * 1000;
}

class CrashAudioProcessor extends AudioWorkletProcessor {
  constructor(options) {
    super();
    this.failed = false;
    this.engine = null;
    this.pendingMessages = [];
    this.port.onmessage = event => {
      if (this.engine) this.#onMessage(event.data);
      else this.pendingMessages.push(event.data);
    };

    const initialMacros = options.processorOptions?.macros ?? [];
    modulePromise.then(module => {
      if (this.failed) return;
      const started = nowMs();
      this.engine = new WasmCrashEngine(module, sampleRate, false);
      this.engine.setMacros(initialMacros);
      this.port.postMessage({
        kind: "ready", elapsedMs: nowMs() - started,
      });
      const pending = this.pendingMessages;
      this.pendingMessages = [];
      pending.forEach(data => this.#onMessage(data));
    }).catch(error => this.#fail(error));
  }

  #onMessage(data) {
    if (this.failed) return;
    try {
      if (data.kind === "macros") {
        const started = nowMs();
        this.engine.setMacros(data.values);
        this.port.postMessage({
          kind: "macros-applied", generation: data.generation,
          elapsedMs: nowMs() - started,
        });
      } else if (data.kind === "trigger") {
        this.engine.trigger(data.event);
        this.port.postMessage({ kind: "triggered" });
      } else if (data.kind === "reset") {
        this.engine.reset();
      }
    } catch (error) {
      this.#fail(error);
    }
  }

  #fail(error) {
    this.failed = true;
    this.port.postMessage({ kind: "error", message: String(error) });
  }

  process(_inputs, outputs) {
    const output = outputs[0]?.[0];
    if (!output) return true;
    if (this.failed || !this.engine) {
      output.fill(0);
      return true;
    }
    try {
      this.engine.processQuantum(output);
    } catch (error) {
      output.fill(0);
      this.#fail(error);
    }
    return true;
  }
}

registerProcessor("triggerfish-crash-renderer", CrashAudioProcessor);
