// Compare an archived pre-change WASM module against the current build.
import {readFile, writeFile} from 'node:fs/promises';
import {pathToFileURL} from 'node:url';
import {resolve} from 'node:path';
import {WasmPercussionEngine} from '../workbench/web/wasm_engine_core.mjs';
import {validateFit, fitMacroValues} from '../workbench/web/state.mjs';
import {recipeAdapter} from '../workbench/web/recipe_adapter.mjs';

const baseline = process.argv[2] ?? 'build/motion-optimization-baseline/triggerfish-percussion.mjs';
const current = 'build/workbench-wasm/triggerfish-percussion.mjs';
const factories = await Promise.all([baseline,current].map(async file =>
  (await import(pathToFileURL(resolve(file)))).default));
const engines = await Promise.all(factories.map(async create => new WasmPercussionEngine(await create(),44100)));
const render = (engine, fit, rate, strength, seconds, hits = 1) => {
  engine.setSampleRate(rate);
  const parameters = fitMacroValues(fit, engine.parameters);
  const routing = recipeAdapter(fit.instrument.recipe).routing(fit.instrument);
  const samples = new Float32Array(Math.round(seconds*rate));
  const start = performance.now();
  engine.setConfiguration(parameters,routing);
  let cursor = 0;
  for (let hit=0;hit<hits;++hit) {
    const frame = Math.round(.5*hit*rate);
    if (frame>cursor) engine.processTo(samples,cursor,frame-cursor);
    engine.trigger({...fit.controls.event,strength,seed:fit.controls.event.seed+hit});
    cursor=frame;
  }
  engine.processTo(samples,cursor,samples.length-cursor);
  return {samples,ms:performance.now()-start};
};
const rows = [];
try {
  for (const instrument of ['gong','crash','ride']) {
    const fit = validateFit(JSON.parse(await readFile(`workbench/web/${instrument}_calibration.fit.json`,'utf8')),engines[0].parameters);
    for (const rate of [44100,48000,96000]) {
      for (const strength of [.3,1]) {
        const pair = engines.map(engine => render(engine,fit,rate,strength,2));
        let maxError = 0, errorEnergy = 0, energy = 0;
        for (let i=0;i<pair[0].samples.length;++i) {
          const x=pair[0].samples[i], y=pair[1].samples[i];
          if (!Number.isFinite(y)) throw Error('Nonfinite optimised audio');
          maxError=Math.max(maxError,Math.abs(x-y)); errorEnergy+=(x-y)**2; energy+=x*x;
        }
        const row = {instrument,rate,strength,baselineMs:pair[0].ms,currentMs:pair[1].ms,
          maxError,errorDb:10*Math.log10(Math.max(1e-30,errorEnergy/Math.max(1e-30,energy)))};
        rows.push(row); console.log(JSON.stringify(row));
        if (maxError>1e-5 || errorEnergy>energy*1e-8) throw Error('Audio equivalence failed');
      }
    }
  }
  const gong = validateFit(JSON.parse(await readFile('workbench/web/gong_calibration.fit.json','utf8')),engines[0].parameters);
  const longPair = engines.map(engine => render(engine,gong,48000,1,8,4));
  const maxError = longPair[0].samples.reduce((error,x,i) => Math.max(error,Math.abs(x-longPair[1].samples[i])),0);
  if (!Number.isFinite(maxError) || maxError>1e-5) throw Error('Eight-second repeated-hit equivalence failed');
  const longTail = {seconds:8,hits:4,rate:48000,maxError};
  console.log(JSON.stringify(longTail));
  await writeFile('build/modal-optimization-comparison.json',JSON.stringify({baseline,current,rows,longTail},null,2));
} finally { for (const engine of engines) engine.destroy(); }
