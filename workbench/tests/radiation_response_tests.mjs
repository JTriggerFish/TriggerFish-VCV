import assert from 'node:assert/strict';
import {radiationCoefficients,responseDb,radiationCurves,spectrumHistogram,eqResponseSampleRate} from '../web/radiation_response.mjs';
const state={reference:{sampleRate:44100},synthesisSpectrum:{sampleRate:44100},
  liveEqHistogram:[-60],liveEqSampleRate:48000};
assert.equal(eqResponseSampleRate(state),48000);
state.liveEqHistogram=null;
assert.equal(eqResponseSampleRate(state),44100);
for(const rate of [44100,48000,96000]) {
  for(const type of ['highpass','lowpass']) {
    const c=radiationCoefficients(type,1000,0,rate);
    assert.ok(Math.abs(responseDb(c,1000,rate)+3.01029995664)<1e-8);
  }
  const c=radiationCoefficients('peak',2500,6,rate);
  assert.ok(Math.abs(responseDb(c,2500,rate)-6)<1e-8);
  const flat=radiationCoefficients('peak',2500,0,rate);
  for(const f of [20,100,1000,10000])assert.ok(Math.abs(responseDb(flat,f,rate))<1e-8);
  const bypass=radiationCurves({low:100,high:15000,frequency:2500,gain:6,enabled:false},rate,[20,1000,19000]);
  assert.deepEqual(bypass.total,[0,0,0]);
}
assert.equal(spectrumHistogram(null),null);
const spectrum={values:new Float32Array(513*20).fill(-60),frames:20,bins:513,size:1024,hop:256,sampleRate:48000};
for(const db of spectrumHistogram(spectrum))assert.ok(db===-160||Math.abs(db-(-60-10*Math.log10(2*1.5*48000/1024)))<1e-8);
console.log('radiation graph: Butterworth cutoffs, peak gain, bypass, fixed spectrum scale');
