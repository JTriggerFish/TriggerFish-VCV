import createModule from "./triggerfish-percussion-audio.mjs";
import { WasmPercussionEngine } from "./wasm_engine_core.mjs";

const modulePromise = Promise.resolve(createModule());

function nowMs() {
  return globalThis.performance?.now() ?? Date.now();
}

class PercussionAudioProcessor extends AudioWorkletProcessor {
  constructor(options) {
    super();
    this.failed = false;
    this.engine = null;
    this.pendingMessages = [];
    this.port.onmessage = event => this.#onMessage(event.data);

    const prepared = options.processorOptions?.prepared;
    const recipeIndex = options.processorOptions?.recipeIndex ?? 0;
    modulePromise.then(module => {
      if (this.failed) return;
      const started = nowMs();
      this.engine = new WasmPercussionEngine(
        module, sampleRate, recipeIndex, false, true,
      );
      this.engine.applyPrepared(prepared);
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
    if (!this.engine) {
      this.pendingMessages.push(data);
      return;
    }
    try {
      if (data.kind === "trigger") {
        this.engine.trigger(data.event);
        this.port.postMessage({ kind: "triggered" });
      } else if (data.kind === "reset") {
        this.engine.reset();
      } else if (data.kind === "dispose") {
        this.engine.destroy();
        this.engine = null;
        this.disposed = true;
        this.port.postMessage({ kind: "disposed" });
      }
    } catch (error) {
      this.#fail(error);
    }
  }

  #fail(error) {
    if (this.failed) return;
    this.failed = true;
    this.engine?.destroy();
    this.engine = null;
    this.pendingMessages = [];
    this.disposed = true;
    this.port.postMessage({ kind: "error", message: String(error) });
  }

  process(_inputs, outputs) {
    const output = outputs[0]?.[0];
    if (!output) return true;
    if (this.failed || !this.engine) {
      output.fill(0);
      return !this.disposed;
    }
    try {
      this.engine.processQuantum(output);
    } catch (error) {
      output.fill(0);
      this.#fail(error);
    }
    return !this.disposed;
  }
}

registerProcessor("triggerfish-percussion-renderer", PercussionAudioProcessor);
