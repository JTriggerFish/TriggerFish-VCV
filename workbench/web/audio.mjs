import { limiterLookaheadSeconds } from "./limiter_config.mjs";

export class SafeAudition {
  constructor(onStatus = () => {}) {
    this.onStatus = onStatus;
    this.context = null;
    this.source = null;
    this.masterDb = -12;
    this.trimDb = 0;
    this.reductionDb = 0;
    this.triggerCount = 0;
    this.rendererReady = false;
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
    this.context = new AudioContext({ latencyHint: "interactive" });
    if (this.context.state === "running") await this.context.suspend();
    await Promise.all([
      this.context.audioWorklet.addModule("crash_audio_worklet_processor.mjs"),
      this.context.audioWorklet.addModule("lookahead_limiter_processor.mjs"),
    ]);

    const initialGeneration = this.macroGeneration;
    this.live = new AudioWorkletNode(
      this.context, "triggerfish-crash-renderer", {
        numberOfInputs: 0, numberOfOutputs: 1, outputChannelCount: [1],
        processorOptions: { macros: [...this.pendingMacros] },
      },
    );
    const ready = this.#bindRenderer(initialGeneration);
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
    await ready;
  }

  #bindRenderer(initialGeneration) {
    return new Promise((resolve, reject) => {
      let starting = true;
      this.live.onprocessorerror = () => {
        const error = new Error("crash AudioWorklet processor failed");
        this.onStatus(error.message);
        if (starting) reject(error);
      };
      this.live.port.onmessage = ({ data }) => {
        if (data.kind === "ready") {
          this.rendererReady = true;
          this.appliedMacroGeneration = initialGeneration;
          this.macroCommitMs = data.elapsedMs;
          starting = false;
          if (this.macroGeneration > initialGeneration) {
            this.live.port.postMessage({
              kind: "macros", values: this.pendingMacros,
              generation: this.macroGeneration,
            });
          }
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
          if (starting) reject(error);
        }
      };
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
    const generation = ++this.macroGeneration;
    this.macroCommitStarted = performance.now();
    if (this.rendererReady) {
      this.live.port.postMessage({ kind: "macros", values, generation });
    }
  }

  async trigger(event) {
    await this.#prepare();
    this.#stopSource();
    this.liveGain.gain.value = 1;
    this.live.port.postMessage({ kind: "trigger", event });
  }

  async activate() { await this.#prepare(); }

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
    this.live?.port.postMessage({ kind: "reset" });
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
  get underflows() { return 0; }
  get outputDb() {
    if (!this.meter) return -Infinity;
    this.meter.getFloatTimeDomainData(this.meterSamples);
    let peak = 0;
    for (const sample of this.meterSamples) peak = Math.max(peak, Math.abs(sample));
    return 20 * Math.log10(Math.max(peak, 1e-9));
  }
  get latencyMs() { return this.latencyBreakdown.totalMs; }
  get latencyBreakdown() {
    const limiterMs = 1000 * limiterLookaheadSeconds;
    if (!this.context) {
      return {
        workletQuantumFrames: 128, workletQuantumMs: 0,
        limiterMs, browserDeviceMs: 0, totalMs: 0,
      };
    }
    const workletQuantumFrames = this.context.renderQuantumSize ?? 128;
    const workletQuantumMs =
      1000 * workletQuantumFrames / this.context.sampleRate;
    const browserDeviceMs = 1000 * ((this.context.baseLatency ?? 0) +
      (this.context.outputLatency ?? 0));
    return {
      workletQuantumFrames, workletQuantumMs, limiterMs, browserDeviceMs,
      totalMs: workletQuantumMs + limiterMs + browserDeviceMs,
    };
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
