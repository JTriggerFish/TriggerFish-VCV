// Actual Wasm integration: no surrogate in acceptance measurements.
import {readFile,writeFile} from 'node:fs/promises';
import {PercussionEngine} from '../build/workbench-wasm/site/engine.mjs';
import {recipeAdapter} from '../workbench/web/recipe_adapter.mjs';
import {decayEnvelope} from '../workbench/web/decay_hold_measurement.mjs';
import {compensateDecay} from '../workbench/web/decay_hold_fit.mjs';
const name=process.argv[2]??'crash';
const fit=JSON.parse(await readFile(`workbench/web/${name}_calibration.fit.json`,'utf8'));
const engine=await PercussionEngine.create(fit.reference.sampleRate,0);
try {
  const raw=Object.assign({},...fit.instrument.nodes.map(n=>n.parameters));
  const baseline=engine.parameters.map(d=>raw[d.key]);
  const edited=baseline.slice(), rate=engine.parameters.find(d=>d.key==='bloom_rate').index;
  edited[rate]*=1.25;
  const start=performance.now();
  const result=await compensateDecay({descriptors:engine.parameters,baseline,edited,seed:fit.controls.event.seed,
    measure:(parameters,seed)=>decayEnvelope(engine.render({seconds:6,parameters,
      routing:recipeAdapter(fit.instrument.recipe).routing(fit.instrument),...fit.controls.event,seed}),engine.sampleRate)});
  const changes=engine.parameters.filter(d=>result.values?.[d.index]!==edited[d.index])
    .map(d=>({key:d.key,before:edited[d.index],after:result.values?.[d.index]}));
  const report={...result,values:undefined,changes,seconds:(performance.now()-start)/1000};
  await writeFile(`build/${name}-hold-decay-audit.json`,JSON.stringify(report,null,2));
  console.log(JSON.stringify(report));
} finally {engine.destroy();}
