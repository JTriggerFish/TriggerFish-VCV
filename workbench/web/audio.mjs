import { limiterLookaheadSeconds } from "./limiter_config.mjs";
import { ConfigurationPreparer } from "./configuration_preparer.mjs";
import { StandbyRenderer } from "./standby_renderer.mjs";

export class SafeAudition {
  constructor(onStatus = () => {}) {
    this.onStatus = onStatus;
    this.masterDb = -12;
    this.trimDb = 0;
    this.reductionDb = 0;
    this.triggerCount = 0;
    this.pendingMacros = [];
    this.pendingRouting = [];
    this.pendingRecipeIndex = 0;
    this.configurationGeneration = 0;
    this.appliedConfigurationGeneration = -1;
    this.macroCommitStarted = 0;
    this.macroCommitMs = 0;
    this.preparationMs = 0;
    this.installationMs = 0;
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
      this.context.audioWorklet.addModule(
        "percussion_audio_worklet_processor.mjs",
      ),
      this.context.audioWorklet.addModule("lookahead_limiter_processor.mjs"),
    ]);
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
    this.trim.connect(this.limiter);
    this.limiter.connect(this.master).connect(this.meter).connect(
      this.context.destination,
    );
    this.renderer = new StandbyRenderer({
      context: this.context,
      destination: this.trim,
      onApplied: configuration => this.#configurationApplied(configuration),
      onTriggered: () => { ++this.triggerCount; },
      onError: error => this.onStatus(String(error)),
    });
    this.preparer = new ConfigurationPreparer({
      onPrepared: configuration => this.renderer.configure(configuration),
      onError: error => this.onStatus(String(error)),
    });
    this.setMaster(this.masterDb);
    this.setTrim(this.trimDb);
    this.#configureRenderer();
    await this.preparer.ready();
    await this.renderer.ready();
  }

  #configurationApplied(configuration) {
    this.appliedConfigurationGeneration = Math.max(
      this.appliedConfigurationGeneration, configuration.generation,
    );
    this.macroCommitMs = configuration.elapsedMs;
    this.preparationMs = configuration.preparationMs;
    this.installationMs = configuration.installationMs;
  }

  #configureRenderer() {
    if (!this.renderer) return;
    this.preparer.configure({
      recipeIndex: this.pendingRecipeIndex,
      values: this.pendingMacros,
      routing: this.pendingRouting,
      sampleRate: this.context.sampleRate,
      generation: this.configurationGeneration,
    });
  }

  #configurationChanged() {
    ++this.configurationGeneration;
    this.macroCommitStarted = performance.now();
    this.#configureRenderer();
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
    this.#configurationChanged();
  }

  setRouting(values) {
    this.pendingRouting = [...values];
    this.#configurationChanged();
  }

  setRecipe(recipeIndex, values, routing) {
    this.pendingRecipeIndex = recipeIndex;
    this.pendingMacros = [...values];
    this.pendingRouting = [...routing];
    this.#configurationChanged();
  }

  async trigger(event) {
    await this.#prepare();
    this.#stopSource();
    await this.preparer.ready();
    await this.renderer.trigger(event);
  }

  async activate() { await this.#prepare(); }

  async play(samples, sampleRate) {
    await this.#prepare();
    this.#stopSource();
    this.renderer.setMuted(true);
    const buffer = new AudioBuffer({
      length: samples.length, numberOfChannels: 1, sampleRate,
    });
    buffer.copyToChannel(samples, 0);
    this.source = new AudioBufferSourceNode(this.context, { buffer });
    this.source.connect(this.trim);
    this.source.onended = () => {
      this.source = null;
      this.renderer.setMuted(false);
    };
    this.source.start();
  }

  stop() {
    this.#stopSource();
    this.renderer?.setMuted(false);
    this.renderer?.reset();
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
    for (const sample of this.meterSamples)
      peak = Math.max(peak, Math.abs(sample));
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
    return this.appliedConfigurationGeneration < this.configurationGeneration;
  }
  get macroCommitElapsedMs() {
    return this.macroCommitPending
      ? performance.now() - this.macroCommitStarted : this.macroCommitMs;
  }
}
