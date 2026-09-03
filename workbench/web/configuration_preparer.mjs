function cloneConfiguration(configuration) {
  return {
    ...configuration,
    values: [...configuration.values],
    routing: [...configuration.routing],
  };
}

export class ConfigurationPreparer {
  constructor({ onPrepared, onError }) {
    this.onPrepared = onPrepared;
    this.onError = onError;
    this.requested = null;
    this.preparedGeneration = -1;
    this.#createWorker();
  }

  #createWorker() {
    this.worker = new Worker("preparation_worker.mjs", { type: "module" });
    this.worker.onerror = event => {
      event.preventDefault();
      const error = new Error(
        event.message || "percussion preparation worker failed",
      );
      const reject = this.rejectCurrent;
      this.rejectCurrent = null;
      this.worker.terminate();
      this.#createWorker();
      if (reject) reject(error); else this.onError(error);
    };
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
    if (!this.requested) return;
    while (this.preparedGeneration < this.requested.generation) {
      if (!this.building) this.configure(this.requested);
      await this.building;
    }
  }

  dispose() {
    this.rejectCurrent?.(new Error("percussion preparation was cancelled"));
    this.rejectCurrent = null;
    this.worker.terminate();
  }

  async #buildLoop() {
    do {
      const configuration = cloneConfiguration(this.requested);
      const result = await this.#prepare(configuration);
      if (configuration.generation !== this.requested.generation) continue;
      this.preparedGeneration = configuration.generation;
      this.onPrepared({
        ...configuration,
        prepared: result.prepared,
        preparationMs: result.elapsedMs,
      });
    } while (this.preparedGeneration < this.requested.generation);
  }

  #prepare(configuration) {
    return new Promise((resolve, reject) => {
      this.worker.onmessage = ({ data }) => {
        if (data.generation !== configuration.generation) return;
        this.rejectCurrent = null;
        if (data.error) reject(new Error(data.error)); else resolve(data);
      };
      this.rejectCurrent = reject;
      this.worker.postMessage(configuration);
    });
  }
}
