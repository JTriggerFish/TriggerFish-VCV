const ringFrames = 8192;
const liveLeadFrames = 512;

export class SafeAudition {
  constructor(onStatus = () => {}) {
    this.onStatus = onStatus;
    this.context = null;
    this.source = null;
    this.masterDb = -12;
    this.trimDb = 0;
    this.reductionDb = 0;
    this.triggerCount = 0;
    this.workerReady = false;
    this.pendingMacros = [];
    this.macroGeneration = 0;
    this.appliedMacroGeneration = 0;
    this.macroCommitStarted = 0;
    this.macroCommitMs = 0;
  }

  async #prepare() {
    await this.initialize();
    await this.context.resume();
  }

  initialize() {
    if (!this.initializing) this.initializing = this.#buildGraph();
    return this.initializing;
  }

  async #buildGraph() {
    if (!crossOriginIsolated || typeof SharedArrayBuffer === "undefined") {
      throw new Error("real-time audio requires the workbench server");
    }
    this.context = new AudioContext({ latencyHint: "interactive" });
    if (this.context.state === "running") await this.context.suspend();
    await Promise.all([
      this.context.audioWorklet.addModule("ring_player_processor.mjs"),
      this.context.audioWorklet.addModule("lookahead_limiter_processor.mjs"),
    ]);
    this.audioBuffer = new SharedArrayBuffer(
      ringFrames * Float32Array.BYTES_PER_ELEMENT,
    );
    this.controlBuffer = new SharedArrayBuffer(4 * Int32Array.BYTES_PER_ELEMENT);
    this.control = new Int32Array(this.controlBuffer);
    this.live = new AudioWorkletNode(this.context, "triggerfish-ring-player", {
      numberOfInputs: 0, numberOfOutputs: 1, outputChannelCount: [1],
      processorOptions: {
        audioBuffer: this.audioBuffer, controlBuffer: this.controlBuffer,
      },
    });
    this.liveGain = new GainNode(this.context);
    this.trim = new GainNode(this.context);
    this.limiter = new AudioWorkletNode(
      this.context, "triggerfish-lookahead-limiter",
      { numberOfInputs: 1, numberOfOutputs: 1, outputChannelCount: [1] },
    );
    this.limiter.port.onmessage = event => {
      this.reductionDb = Number(event.data);
    };
    this.master = new GainNode(this.context);
    this.meter = new AnalyserNode(this.context, {
      fftSize: 256, smoothingTimeConstant: 0,
    });
    this.meterSamples = new Float32Array(this.meter.fftSize);
    this.live.connect(this.liveGain).connect(this.trim).connect(this.limiter);
    this.limiter.connect(this.master).connect(this.meter).connect(
      this.context.destination,
    );
    this.setMaster(this.masterDb);
    this.setTrim(this.trimDb);
    await this.#startWorker();
    Atomics.store(this.control, 2, 0);
  }

  #startWorker() {
    this.worker = new Worker("realtime_engine_worker.mjs", { type: "module" });
    return new Promise((resolve, reject) => {
      this.worker.onmessage = ({ data }) => {
        if (data.kind === "ready") {
          this.workerReady = true;
          const generation = ++this.macroGeneration;
          this.macroCommitStarted = performance.now();
          this.worker.postMessage({
            kind: "macros", values: this.pendingMacros, generation,
          });
          resolve();
        } else if (data.kind === "macros-applied") {
          if (data.generation >= this.appliedMacroGeneration) {
            this.appliedMacroGeneration = data.generation;
            this.macroCommitMs = data.elapsedMs;
          }
        } else if (data.kind === "triggered") {
          ++this.triggerCount;
        } else if (data.kind === "error") {
          const error = new Error(data.message);
          this.onStatus(error.message);
          reject(error);
        }
      };
      this.worker.onerror = event => reject(new Error(event.message));
      this.worker.postMessage({
        kind: "initialize", sampleRate: this.context.sampleRate,
        macros: this.pendingMacros, audioBuffer: this.audioBuffer,
        controlBuffer: this.controlBuffer,
      });
    });
  }

  setMaster(db) {
    this.masterDb = Math.min(0, Math.max(-60, Number(db)));
    if (this.master) this.master.gain.value = 10 ** (this.masterDb / 20);
  }

  setTrim(db) {
    this.trimDb = Math.min(36, Math.max(-36, Number(db)));
    if (this.trim) this.trim.gain.value = 10 ** (this.trimDb / 20);
  }

  setMacros(values) {
    this.pendingMacros = [...values];
    if (this.workerReady) {
      const generation = ++this.macroGeneration;
      this.macroCommitStarted = performance.now();
      this.worker.postMessage({ kind: "macros", values, generation });
    }
  }

  async trigger(event) {
    await this.#prepare();
    this.#stopSource();
    this.liveGain.gain.value = 1;
    this.worker.postMessage({ kind: "trigger", event });
  }

  async play(samples, sampleRate) {
    await this.#prepare();
    this.#stopSource();
    this.liveGain.gain.value = 0;
    const buffer = new AudioBuffer({
      length: samples.length, numberOfChannels: 1, sampleRate,
    });
    buffer.copyToChannel(samples, 0);
    this.source = new AudioBufferSourceNode(this.context, { buffer });
    this.source.connect(this.trim);
    this.source.onended = () => {
      this.source = null;
      this.liveGain.gain.value = 1;
    };
    this.source.start();
  }

  stop() {
    this.#stopSource();
    if (this.liveGain) this.liveGain.gain.value = 1;
    this.worker?.postMessage({ kind: "reset" });
  }

  #stopSource() {
    if (!this.source) return;
    this.source.onended = null;
    try { this.source.stop(); } catch { /* already stopped */ }
    this.source.disconnect();
    this.source = null;
  }

  get reduction() { return this.reductionDb; }
  get state() { return this.context?.state ?? "off"; }
  get underflows() { return this.control ? Atomics.load(this.control, 2) : 0; }
  get outputDb() {
    if (!this.meter) return -Infinity;
    this.meter.getFloatTimeDomainData(this.meterSamples);
    let peak = 0;
    for (const sample of this.meterSamples) peak = Math.max(peak, Math.abs(sample));
    return 20 * Math.log10(Math.max(peak, 1e-9));
  }
  get latencyMs() {
    if (!this.context) return 0;
    return 1000 * ((this.context.baseLatency ?? 0) +
      (this.context.outputLatency ?? 0) +
      liveLeadFrames / this.context.sampleRate + 0.005);
  }
  get sampleRate() { return this.context?.sampleRate ?? 0; }
  get macroCommitPending() {
    return this.appliedMacroGeneration < this.macroGeneration;
  }
  get macroCommitElapsedMs() {
    return this.macroCommitPending
      ? performance.now() - this.macroCommitStarted : this.macroCommitMs;
  }
}
