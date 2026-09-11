// Silent exact-WASM snapshot audit through the actual browser safety limiter.
import { readFile, mkdir, writeFile } from 'node:fs/promises';
import { resolve } from 'node:path';
import { pathToFileURL } from 'node:url';

const [filename, destination = 'build/limiter-fit-audit', trialIndex] = process.argv.slice(2);
if (!filename) throw Error('Usage: node tools/audit_fit_limiter.mjs fit.json [output-directory] [trial-index]');
const load = name => import(pathToFileURL(resolve('build/workbench-wasm/site', name)));
const { PercussionEngine } = await load('engine.mjs');
const { validateFit, fitMacroValues } = await load('state.mjs');
const { recipeAdapter } = await load('recipe_adapter.mjs');
globalThis.sampleRate = 48000;
globalThis.AudioWorkletProcessor = class {
  constructor() { this.port = { postMessage() {} }; }
};
let Limiter;
globalThis.registerProcessor = (_name, type) => { Limiter = type; };
await import('../workbench/web/lookahead_limiter_processor.mjs');

function measure(samples, masterDb) {
  const limiter = new Limiter(), gain = 10 ** (masterDb / 20);
  let rawPeak = 0, outputPeak = 0, reduction = 0, limitedFrames = 0;
  for (let first = 0; first < samples.length; first += 128) {
    const input = Float32Array.from(samples.subarray(first, first + 128), x => x * gain);
    const output = new Float32Array(input.length);
    limiter.process([[input]], [[output]]);
    reduction = Math.min(reduction, limiter.gainDb);
    if (limiter.gainDb < -.1) limitedFrames += input.length;
    for (const x of input) rawPeak = Math.max(rawPeak, Math.abs(x));
    for (const x of output) outputPeak = Math.max(outputPeak, Math.abs(x));
  }
  return { masterDb, inputSamplePeakDb: 20 * Math.log10(rawPeak),
    outputSamplePeakDb: 20 * Math.log10(outputPeak), maximumReductionDb: -reduction,
    limitedSecondsApprox: limitedFrames / sampleRate };
}

const engine = await PercussionEngine.create(sampleRate);
try {
  const document = JSON.parse(await readFile(filename, 'utf8'));
  const index = trialIndex === undefined ? undefined : Number(trialIndex);
  if (Array.isArray(document) && (!Number.isInteger(index) || index < 0 || index >= document.length))
    throw Error('A fit array requires a valid zero-based trial index');
  if (!Array.isArray(document) && index !== undefined)
    throw Error('A trial index is only valid for a fit array');
  const source = Array.isArray(document) ? document[index] : document;
  engine.setRecipe(engine.recipes.find(x => x.key === source.instrument.recipe).index);
  const fit = validateFit(source, engine.parameters);
  const values = fitMacroValues(fit, engine.parameters);
  const routing = recipeAdapter(fit.instrument.recipe).routing(fit.instrument);
  const rows = [];
  for (const hits of [1, 4]) {
    engine.reset();
    engine.setConfiguration(values, routing);
    const samples = new Float32Array(sampleRate * 6);
    let cursor = 0;
    for (let hit = 0; hit < hits; ++hit) {
      const frame = Math.round(hit * .5 * sampleRate);
      if (frame > cursor) engine.processTo(samples, cursor, frame - cursor);
      engine.trigger({ ...fit.controls.event, seed: (fit.controls.event.seed + hit) >>> 0 });
      cursor = frame;
    }
    engine.processTo(samples, cursor, samples.length - cursor);
    if (!samples.every(Number.isFinite)) throw Error('Nonfinite synthesis');
    rows.push({ hits, measurements: [-12, 0].map(db => measure(samples, db)) });
  }
  const report = { source: resolve(filename), trialIndex: index, name: fit.name, sampleRate,
    event: fit.controls.event, notes: 'Current UI snapshot migration; actual limiter, no playback. Master is not saved in fits.', rows };
  await mkdir(destination, { recursive: true });
  await writeFile(resolve(destination, 'audit.json'), JSON.stringify(report, null, 2));
  console.log(JSON.stringify(report, null, 2));
} finally { engine.destroy(); }
