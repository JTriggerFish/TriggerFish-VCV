import assert from "node:assert/strict";
import { LimiterMeter, paintLimiterMeter } from '../web/limiter_meter.mjs';

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
assert.equal(limiter.lookahead, 144);
const input = new Float32Array(4096);
input[64] = 2;
const output = process(limiter, input);
const peak = Math.max(...output.map(Math.abs));
const peakIndex = output.findIndex(value => Math.abs(value) === peak);
assert.ok(peak <= 10 ** (-1 / 20) + 1e-6, `unsafe peak ${peak}`);
assert.equal(peakIndex, 64 + limiter.lookahead);

// Report the interval maximum, not only gain at the last reporting sample.
const metered = new Limiter(), messages = [];
metered.port.postMessage = value => messages.push(value);
process(metered, input);
assert.equal(messages.length, 2);
assert.ok(messages[0].intervalReductionDb < -6);
assert.ok(messages[0].intervalReductionDb <= messages[0].reductionDb);
assert.ok(messages[0].inputPeakDb > 6);

const display = new LimiterMeter();
display.ingest({ reductionDb: 0, intervalReductionDb: -.5 }, 100);
assert.deepEqual(display.read(200), {active: true, recentDb: .5, maximumDb: .5});
display.ingest({ reductionDb: 0, intervalReductionDb: 0 }, 300);
assert.equal(display.read(1050).active, true);
assert.deepEqual(display.read(1100), {active: false, recentDb: 0, maximumDb: .5});
display.ingest({ reductionDb: -4 }, 1200);
display.ingest({ reductionDb: NaN }, 1300);
assert.equal(display.read(1300).maximumDb, 4);
display.reset();
assert.deepEqual(display.read(1400), {active: false, recentDb: 0, maximumDb: 0});

const classes = new Set(), nodes = new Map();
for (const key of ['.limiter-label', 'meter', 'output', 'button'])
  nodes.set(key, {max: 24, setAttribute() {}});
const element = { classList: {toggle(key, on) {on ? classes.add(key) : classes.delete(key);}},
  querySelector(key) {return nodes.get(key);} };
paintLimiterMeter(element, {active: true, recentDb: 10.7, maximumDb: 10.7});
assert.equal(nodes.get('.limiter-label').textContent, 'LIMITING');
assert.equal(nodes.get('meter').value, 10.7);
assert.ok(classes.has('is-limiting'));
paintLimiterMeter(element, {active: false, recentDb: 0, maximumDb: 10.7});
assert.ok(!classes.has('is-limiting') && classes.has('has-limited'));
assert.equal(nodes.get('button').textContent, 'Max 10.7 dB · reset');

const unityLimiter = new Limiter();
const quiet = new Float32Array(2048);
quiet[32] = 0.25;
const quietOutput = process(unityLimiter, quiet);
assert.ok(Math.abs(quietOutput[32 + unityLimiter.lookahead] - 0.25) < 1e-7);

function peakOf(samples) {
  return samples.reduce((peak, value) => Math.max(peak, Math.abs(value)), 0);
}

// Independent longer-kernel, 8x reconstruction, not the production detector.
function reconstructedPeak(samples) {
  const radius = 32;
  let peak = peakOf(samples);
  for (let phase = 1; phase < 8; ++phase) {
    const kernel = Float64Array.from({ length: 2*radius + 1 }, (_, tap) => {
      const x = tap - radius + phase/8;
      const window = Math.abs(x) < radius
        ? .42 + .5*Math.cos(Math.PI*x/radius) + .08*Math.cos(2*Math.PI*x/radius) : 0;
      return window * Math.sin(Math.PI*x) / (Math.PI*x);
    });
    const norm = kernel.reduce((a, b) => a+b, 0);
    for (let sample = radius; sample < samples.length - radius; ++sample) {
      let value = 0;
      for (let tap = 0; tap < kernel.length; ++tap)
        value += kernel[tap] * samples[sample + radius - tap];
      peak = Math.max(peak, Math.abs(value / norm));
    }
  }
  return peak;
}

for (const rate of [44100, 48000, 96000]) {
  globalThis.sampleRate = rate;
  for (const frequency of [20, 40, 80, 128, 250, rate/4, rate*.43]) {
    for (const amplitude of [.1, 4, 350, 1e6]) {
      const processor = new Limiter();
      const signal = Float32Array.from({ length: Math.ceil(.18*rate) },
        (_, sample) => amplitude * Math.sin(2*Math.PI*frequency*sample/rate + .7));
      const rendered = process(processor, signal);
      assert.ok(peakOf(rendered) <= 10**(-1/20) + 1e-6,
        `unsafe sustained sample peak at ${frequency} Hz, A=${amplitude}`);
      if (amplitude === .1) {
        for (let sample = processor.lookahead; sample < rendered.length; ++sample)
          assert.equal(rendered[sample], signal[sample - processor.lookahead]);
      }
      if (amplitude === 350 && frequency < 1000) {
        let error = 0, energy = 0;
        for (let sample = Math.ceil(.12*rate); sample < rendered.length; ++sample) {
          const expected = signal[sample - processor.lookahead] * processor.gain;
          error += (rendered[sample] - expected)**2;
          energy += expected**2;
        }
        assert.ok(error / energy < 1e-8,
          `limiter distorts steady bass at ${frequency} Hz: ${error/energy}`);
      }
      if (amplitude === 4 && frequency >= rate/4) {
        const reconstructed = reconstructedPeak(rendered);
        assert.ok(reconstructed <= 10**(-1/20),
          `unsafe reconstructed peak ${reconstructed} at ${frequency} Hz`);
      }
    }
  }
}
console.log("lookahead limiter: impulse, sustained overload, bass purity, unity and reconstructed peaks passed");
