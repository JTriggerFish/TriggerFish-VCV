// Diagnostic ablations only. No preset writes or audio playback.
import {readFile} from 'node:fs/promises';
import {PercussionEngine} from '../build/workbench-wasm/site/engine.mjs';
import {fitMacroValues, validateFit} from '../build/workbench-wasm/site/state.mjs';
import {recipeAdapter} from '../build/workbench-wasm/site/recipe_adapter.mjs';

const engine = await PercussionEngine.create(44100);
try {
  const fit = validateFit(JSON.parse(await readFile('workbench/web/gong_calibration.fit.json', 'utf8')), engine.parameters);
  const base = fitMacroValues(fit, engine.parameters);
  const routing = recipeAdapter(fit.instrument.recipe).routing(fit.instrument);
  const cases = [ ['Current gong', {}], ['Movement bypassed', {field_motion_depth: 0}],
    ['Diffusion bypassed', {bloom_rate: 0}], ['Both bypassed', {field_motion_depth: 0, bloom_rate: 0}] ];
  const rows = [];
  for (const [name, changes] of cases) {
    const parameters = [...base];
    for (const [key, value] of Object.entries(changes)) {
      const index = engine.parameters.findIndex(p => p.key === key);
      if (index < 0) throw Error('Unknown parameter: ' + key);
      parameters[index] = value;
    }
    const times = [];
    for (let run = 0; run < 4; ++run) {
      const start = performance.now();
      engine.render({seconds: 2, parameters, routing, ...fit.controls.event});
      if (run) times.push(performance.now() - start);
    }
    const ms = times.sort((a,b) => a-b)[1];
    rows.push({name, medianMs: ms, timesRealTime: 2000 / ms});
  }
  console.log(JSON.stringify(rows, null, 2));
} finally { engine.destroy(); }
