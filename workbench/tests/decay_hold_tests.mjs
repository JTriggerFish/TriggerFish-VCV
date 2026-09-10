import assert from 'node:assert/strict';
import {holdTargets,decayEnvelope} from '../web/decay_hold_measurement.mjs';
import {compensationCoordinates,leastSquaresStep} from '../web/decay_hold_solver.mjs';
import {compensateDecay} from '../web/decay_hold_fit.mjs';

const descriptors = [
  ['body_decay_seconds_0',.1,30], ['body_decay_seconds_7',.1,30],
  ['bloom_rate',0,16], ['bloom_energy_acceleration',0,1], ['model_level_db',-60,12],
].map(([key,minimum,maximum],index)=>({key,minimum,maximum,index}));
const baseline=[4,1,2,.2,-20], edited=[4,1,2.5,.2,-20];
assert.ok(!compensationCoordinates([...descriptors,
  {key:'bloom_energy_sensitivity',index:5,minimum:0,maximum:2}],
  [...baseline,.4],[...edited,.4]).some(x=>x.d.key==='bloom_energy_sensitivity'),
  'Hold decay must not retune velocity sensitivity');
assert.ok(Math.abs(leastSquaresStep([[1,0],[0,1]],[2,3])[0]+2)<.001);
assert.deepEqual(compensationCoordinates(descriptors,baseline,edited).map(a=>a.d.index),[0,1,3]);
assert.deepEqual(compensationCoordinates(descriptors,baseline,[4,1,2.5,.3,-20]).map(a=>a.d.index),[0,1]);

// Known surrogate only for optimizer unit tests, not fitting real instruments.
const measure = async p => Float64Array.from({length:264},(_,i)=> {
  const db=i<120?p[2]: (i%2?6*Math.log(p[0]/4):6*Math.log(p[1]))+4*(p[2]-2);
  return 10**(db/10);
});
const result=await compensateDecay({descriptors,baseline,edited,seed:7,measure});
assert.ok(result.accepted,JSON.stringify(result));
assert.ok(result.after<.1);
assert.equal(result.values[2],edited[2]); // User's transport choice never undone.
assert.equal(result.values[4],-20); // No gain matching.
assert.equal(result.front,0); // Preserve intended new attack, not the old one.
assert.deepEqual(baseline,[4,1,2,.2,-20]);
assert.deepEqual(edited,[4,1,2.5,.2,-20]);
const stable=await compensateDecay({descriptors,baseline,edited:baseline,seed:7,measure});
assert.equal(stable.accepted,false);
const limits=[.1,.1,2,.2,-20];
const rejected=await compensateDecay({descriptors,baseline:limits,edited:[.1,.1,16,.2,-20],seed:7,measure});
assert.equal(rejected.accepted,false);
assert.deepEqual(rejected.values,[.1,.1,16,.2,-20]);
assert.throws(()=>holdTargets(new Float64Array(264),new Float64Array(264)),/audible/);
assert.throws(()=>holdTargets([NaN],[1]),/Invalid/);

const fs=8000;
const tone=Float32Array.from({length:6*fs},(_,i)=>[150,400,1000,2200]
  .reduce((s,f)=>s+Math.sin(2*Math.PI*f*i/fs),0)*Math.exp(-.5*i/fs));
const envelope=decayEnvelope(tone,fs);
assert.equal(envelope.length,264);
assert.ok(envelope.every(Number.isFinite));
const changes=holdTargets(envelope,envelope).errors(envelope);
assert.ok(changes.late.every(x=>x===0));
assert.throws(()=>decayEnvelope(tone.slice(0,100),fs),/six-second/);
console.log('Hold decay: bounded solve, preserved edit/attack/gain, rejection and measurement pass');
