import assert from 'node:assert/strict';
import {stft} from '../web/analysis.mjs';
import {ProgressiveStft} from '../web/progressive_stft.mjs';

const audio = Float32Array.from({length: 8000}, (_, i) =>
  .3 * Math.sin(i * .33) + .1 * Math.cos(i * i * .001));
for (const window of ['hann', 'blackman-harris', 'rectangular']) {
  const rolling = new ProgressiveStft();
  const settings = {size: 256, hop: 64, window};
  const expected = stft(audio, 8000, settings);
  for (const length of [50, 512, 1536, 3000, 6000]) {
    const result = rolling.analyze(audio.slice(0, length), 8000, settings, 1, true);
    assert.deepEqual(result.values, expected.values.slice(0, result.values.length));
    assert.equal(result.incomplete, true);
  }
  const complete = rolling.analyze(audio, 8000, settings, 1, false);
  assert.deepEqual(complete.values, expected.values);
  assert.equal(complete.peakDb, expected.peakDb);
  assert.equal(rolling.cache, null);
  rolling.analyze(audio.slice(0, 2000), 8000, settings, 2, true);
  const silent = rolling.analyze(new Float32Array(8000), 8000, settings, 3, false);
  assert.ok(silent.values.every(x => x === silent.floorDb), 'New render discards old frames');
  rolling.analyze(audio.slice(0, 2000), 8000, settings, 4, true);
  const changed = {...settings, size: 512, hop: 128};
  assert.deepEqual(rolling.analyze(audio, 8000, changed, 4, false).values,
    stft(audio, 8000, changed).values, 'FFT setting changes invalidate cache');
}
console.log('Progressive STFT: exact full-render equivalence, windows and cache invalidation pass');
