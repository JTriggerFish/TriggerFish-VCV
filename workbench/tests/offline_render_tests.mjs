import assert from 'node:assert/strict';
import {pathToFileURL} from 'node:url';
import {resolve} from 'node:path';
import {renderOffline} from '../web/offline_render.mjs';
import {WasmPercussionEngine} from '../web/wasm_engine_core.mjs';

let time = 0, cursor = 0, yields = 0;
const engine = {
  setConfiguration() { cursor = 0; }, trigger() {},
  processTo(array, offset, count) {
    for (let i = 0; i < count; ++i) array[offset+i] = cursor++;
    time += 15;
  },
};
const request = {seconds: 1, sampleRate: 48000, parameters: [], event: {}};
const previews = [];
const result = await renderOffline(engine, request, {
  now: () => time, yieldTask: async () => {++yields;},
  progress: samples => previews.push(samples),
});
assert.equal(result.length, 48000);
assert.ok(yields > 0 && previews.length > 0);
for (const prefix of previews) assert.deepEqual(prefix, result.slice(0, prefix.length));
assert.equal(result[47999], 47999);
let cancel = false;
assert.equal(await renderOffline(engine, request, {
  now: () => time, yieldTask: async () => {cancel = true;},
  cancelled: () => cancel,
}), null);
assert.ok(cursor < 48000, 'Obsolete render stops before its tail');

if (process.argv[2]) {
  const {default: createModule} = await import(pathToFileURL(resolve(process.argv[2])));
  const dsp = new WasmPercussionEngine(await createModule(), 44100, 0);
  try {
    const parameters = dsp.parameters.map(p => p.defaultValue);
    const request = {seconds: .8, sampleRate: 44100, parameters, event: {seed: 1675}};
    dsp.setConfiguration(parameters);
    dsp.trigger(request.event);
    const expected = new Float32Array(Math.round(request.seconds * request.sampleRate));
    dsp.processTo(expected, 0, expected.length);
    const actual = await renderOffline(dsp, request);
    assert.deepEqual(actual, expected, 'Chunking must not change the actual DSP output');
  } finally { dsp.destroy(); }
}
console.log('Offline render: exact prefixes, cancellation and complete output pass');
