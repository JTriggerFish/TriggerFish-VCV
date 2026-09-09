import {PercussionEngine} from './engine.mjs';
import {decayEnvelope} from './decay_hold_measurement.mjs';
import {compensateDecay} from './decay_hold_fit.mjs';

self.onmessage=async({data})=>{
  let engine;
  try {
    const {sampleRate,recipeIndex,routing,event,baseline,edited}=data;
    engine=await PercussionEngine.create(sampleRate,recipeIndex);
    const result=await compensateDecay({baseline,edited,descriptors:engine.parameters,seed:event.seed,
      measure:(parameters,seed)=>decayEnvelope(engine.render({seconds:6,parameters,routing,...event,seed}),sampleRate),
      progress:evaluations=>self.postMessage({progress:evaluations})});
    self.postMessage({result});
  } catch(error) { self.postMessage({error:String(error)}); }
  finally {engine?.destroy();}
};
