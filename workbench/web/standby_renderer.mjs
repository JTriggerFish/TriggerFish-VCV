const CrossfadeSeconds = .003;

function cloneConfiguration(configuration) {
  return {
    ...configuration,
    values: [...configuration.values],
    routing: [...configuration.routing],
    prepared: configuration.prepared?.slice(0),
  };
}

export class StandbyRenderer {
  constructor({ context, destination, onApplied, onTriggered, onError }) {
    this.context = context;
    this.destination = destination;
    this.onApplied = onApplied;
    this.onTriggered = onTriggered;
    this.onError = onError;
    this.requested = null;
    this.active = null;
    this.appliedGeneration = -1;
    this.muted = false;
    this.building = null;
    this.retirement = Promise.resolve();
  }

  configure(configuration) {
    this.requested = cloneConfiguration(configuration);
    if (!this.building) {
      this.building = this.#buildLoop().catch(error => {
        this.onError(error);
        throw error;
      }).finally(() => { this.building = null; });
      this.building.catch(() => {});
    }
  }

  async ready() {
    while (!this.active || this.appliedGeneration < this.requested.generation) {
      if (!this.building) this.configure(this.requested);
      await this.building;
    }
  }

  async trigger(event) {
    await this.ready();
    this.setMuted(false);
    this.active.node.port.postMessage({ kind: "trigger", event });
  }

  reset() {
    this.active?.node.port.postMessage({ kind: "reset" });
  }

  setMuted(muted) {
    this.muted = Boolean(muted);
    if (!this.active) return;
    const gain = this.active.gain.gain;
    gain.cancelScheduledValues(this.context.currentTime);
    gain.setValueAtTime(this.muted ? 0 : 1, this.context.currentTime);
  }

  async #buildLoop() {
    do {
      // Keep at most the active and one candidate Wasm session alive. This
      // also bounds rapid slider gestures when the fixed session pool is full.
      await this.retirement;
      const configuration = cloneConfiguration(this.requested);
      const candidate = await this.#create(configuration);
      if (configuration.generation !== this.requested.generation) {
        await this.#dispose(candidate);
        continue;
      }
      this.#activate(candidate, configuration);
    } while (this.appliedGeneration < this.requested.generation);
  }

  #create(configuration) {
    return new Promise((resolve, reject) => {
      const node = new AudioWorkletNode(
        this.context, "triggerfish-percussion-renderer", {
          numberOfInputs: 0, numberOfOutputs: 1, outputChannelCount: [1],
          processorOptions: configuration,
        },
      );
      const gain = new GainNode(this.context, { gain: 0 });
      const renderer = { node, gain, disposed: false };
      node.connect(gain).connect(this.destination);
      let starting = true;
      const fail = error => {
        if (renderer.disposed) return;
        node.disconnect();
        gain.disconnect();
        node.port.close();
        renderer.disposed = true;
        if (starting) {
          starting = false;
          reject(error);
        } else {
          if (this.active === renderer) this.active = null;
          this.onError(error);
        }
      };
      node.onprocessorerror = () => {
        fail(new Error("percussion AudioWorklet processor failed"));
      };
      node.port.onmessage = ({ data }) => {
        if (data.kind === "ready") {
          starting = false;
          renderer.elapsedMs = data.elapsedMs;
          resolve(renderer);
        } else if (data.kind === "triggered") {
          this.onTriggered();
        } else if (data.kind === "disposed") {
          renderer.finishDisposal?.();
        } else if (data.kind === "error") {
          fail(new Error(data.message));
        }
      };
    });
  }

  #activate(candidate, configuration) {
    const previous = this.active;
    this.active = candidate;
    this.appliedGeneration = configuration.generation;
    const now = this.context.currentTime;
    const target = this.muted ? 0 : 1;
    candidate.gain.gain.setValueAtTime(0, now);
    candidate.gain.gain.linearRampToValueAtTime(
      target, now + CrossfadeSeconds,
    );
    if (previous) {
      previous.gain.gain.cancelScheduledValues(now);
      previous.gain.gain.setValueAtTime(previous.gain.gain.value, now);
      previous.gain.gain.linearRampToValueAtTime(0, now + CrossfadeSeconds);
      this.retirement = new Promise(resolve => {
        setTimeout(() => {
          this.#dispose(previous).finally(resolve);
        }, 1000 * CrossfadeSeconds + 8);
      });
    }
    this.onApplied({
      ...configuration,
      installationMs: candidate.elapsedMs,
      elapsedMs: configuration.preparationMs + candidate.elapsedMs,
    });
  }

  #dispose(renderer) {
    if (renderer.disposal) return renderer.disposal;
    renderer.disposal = new Promise(resolve => {
      const finish = () => {
        if (renderer.disposed) {
          resolve();
          return;
        }
        renderer.disposed = true;
        clearTimeout(timeout);
        renderer.node.disconnect();
        renderer.gain.disconnect();
        renderer.node.port.close();
        resolve();
      };
      renderer.finishDisposal = finish;
      const timeout = setTimeout(finish, 100);
      renderer.node.port.postMessage({ kind: "dispose" });
    });
    return renderer.disposal;
  }
}
