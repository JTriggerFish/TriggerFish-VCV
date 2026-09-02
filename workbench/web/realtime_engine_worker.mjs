import { CrashEngine } from "./engine.mjs";

const ReadIndex = 0;
const WriteIndex = 1;
const Underflows = 2;
const blockFrames = 256;
const targetLeadFrames = 512;
let engine;
let audio;
let control;
let capacity;
let running = false;

function queuedFrames() {
  const read = Atomics.load(control, ReadIndex) >>> 0;
  const write = Atomics.load(control, WriteIndex) >>> 0;
  return (write - read) >>> 0;
}

function pump() {
  if (!engine || queuedFrames() >= targetLeadFrames) return;
  const write = Atomics.load(control, WriteIndex) >>> 0;
  const count = Math.min(blockFrames, targetLeadFrames - queuedFrames());
  const first = write % capacity;
  const beforeWrap = Math.min(count, capacity - first);
  engine.processTo(audio, first, beforeWrap);
  if (beforeWrap < count) engine.processTo(audio, 0, count - beforeWrap);
  Atomics.store(control, WriteIndex, (write + count) | 0);
}

function prime() {
  while (queuedFrames() < targetLeadFrames) pump();
}

async function producerLoop() {
  while (running) {
    const observedRead = Atomics.load(control, ReadIndex);
    prime();
    if (typeof Atomics.waitAsync === "function") {
      await Atomics.waitAsync(control, ReadIndex, observedRead, 100).value;
    } else {
      await new Promise(resolve => setTimeout(resolve, 1));
    }
  }
}

self.onmessage = async ({ data }) => {
  try {
    if (data.kind === "initialize") {
      audio = new Float32Array(data.audioBuffer);
      control = new Int32Array(data.controlBuffer);
      capacity = audio.length;
      engine = await CrashEngine.create(data.sampleRate);
      engine.setMacros(data.macros);
      running = true;
      prime();
      producerLoop();
      self.postMessage({ kind: "ready" });
    } else if (data.kind === "macros") {
      const started = performance.now();
      engine.setMacros(data.values);
      self.postMessage({
        kind: "macros-applied", generation: data.generation,
        elapsedMs: performance.now() - started,
      });
    } else if (data.kind === "trigger") {
      engine.trigger(data.event);
      self.postMessage({ kind: "triggered" });
    } else if (data.kind === "reset") {
      engine.reset();
      Atomics.store(control, Underflows, 0);
    }
  } catch (error) {
    running = false;
    self.postMessage({ kind: "error", message: String(error) });
  }
};
