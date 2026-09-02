import assert from "node:assert/strict";

globalThis.AudioWorkletProcessor = class {};
let RingPlayer;
globalThis.registerProcessor = (_name, processor) => { RingPlayer = processor; };
await import("../web/ring_player_processor.mjs");

const audioBuffer = new SharedArrayBuffer(16 * Float32Array.BYTES_PER_ELEMENT);
const controlBuffer = new SharedArrayBuffer(4 * Int32Array.BYTES_PER_ELEMENT);
const audio = new Float32Array(audioBuffer);
const control = new Int32Array(controlBuffer);
audio.set([0.1, 0.2, 0.3, 0.4]);
control[1] = 4;

const player = new RingPlayer({ processorOptions: { audioBuffer, controlBuffer } });
const output = new Float32Array(6);
assert.equal(player.process([], [[output]]), true);
assert.deepEqual([...output].map(value => value.toFixed(1)),
  ["0.1", "0.2", "0.3", "0.4", "0.0", "0.0"]);
assert.equal(control[0], 4);
assert.equal(control[2], 1);
console.log("ring player tests passed");
