import assert from "node:assert/strict";

globalThis.sampleRate = 48000;
globalThis.AudioWorkletProcessor = class {
  constructor() { this.port = { postMessage() {} }; }
};
let Limiter;
globalThis.registerProcessor = (_name, processor) => { Limiter = processor; };
await import("../web/lookahead_limiter_processor.mjs");

function process(processor, samples) {
  const rendered = [];
  for (let first = 0; first < samples.length; first += 128) {
    const input = samples.slice(first, first + 128);
    const output = new Float32Array(input.length);
    processor.process([[input]], [[output]]);
    rendered.push(...output);
  }
  return rendered;
}

const limiter = new Limiter();
const input = new Float32Array(4096);
input[64] = 2;
const output = process(limiter, input);
const peak = Math.max(...output.map(Math.abs));
const peakIndex = output.findIndex(value => Math.abs(value) === peak);
assert.ok(peak <= 10 ** (-1 / 20) + 1e-6, `unsafe peak ${peak}`);
assert.equal(peakIndex, 64 + limiter.lookahead);

const unityLimiter = new Limiter();
const quiet = new Float32Array(2048);
quiet[32] = 0.25;
const quietOutput = process(unityLimiter, quiet);
assert.ok(Math.abs(quietOutput[32 + unityLimiter.lookahead] - 0.25) < 1e-7);
console.log("lookahead limiter tests passed");
