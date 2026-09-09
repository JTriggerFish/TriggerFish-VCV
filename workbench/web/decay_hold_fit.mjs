import {holdTargets} from './decay_hold_measurement.mjs';
import {compensationCoordinates,leastSquaresStep,rms} from './decay_hold_solver.mjs';

// measure(values, seed) returns power cells from the actual C++ renderer.
export async function compensateDecay({descriptors,baseline,edited,seed,measure,progress=()=>{}}) {
  let evaluations=0;
  const probe=async(values,s=seed)=>{
    progress(++evaluations);
    return await measure(values,s);
  };
  const target=holdTargets(await probe(baseline),await probe(edited));
  const axes=compensationCoordinates(descriptors,baseline,edited);
  const origin=axes.map(a=>a.origin);
  const unpack=x=>{
    const p=edited.slice();axes.forEach((a,i)=>p[a.d.index]=Math.max(a.d.minimum,Math.min(a.d.maximum,a.decode(x[i]))));
    return p;
  };
  const cache=new Map();
  async function evaluate(x) {
    const key=x.join(',');if(cache.has(key))return cache.get(key);
    const errors=target.errors(await probe(unpack(x)));
    const residual=[...errors.late.map(v=>v/Math.sqrt(errors.late.length)),
      ...errors.front.map(v=>3*Math.sign(v)*Math.max(0,Math.abs(v)-1)/Math.sqrt(errors.front.length)),
      ...x.map((v,i)=>.03*(v-origin[i])/axes[i].maxStep)];
    const row={x:x.slice(),errors,residual,score:residual.reduce((s,v)=>s+v*v,0)};
    cache.set(key,row);return row;
  }
  let best=await evaluate(origin);
  const before=rms(best.errors.late);
  if(before<.15)return {accepted:false,reason:'Decay already stable',before,after:before,evaluations};
  for(let iteration=0;iteration<3;iteration++) {
    const columns=[];
    for(let i=0;i<axes.length;i++) {
      const a=best.x.slice(),b=best.x.slice(),axis=axes[i];
      a[i]=Math.max(axis.lo,a[i]-axis.step);b[i]=Math.min(axis.hi,b[i]+axis.step);
      const minus=await evaluate(a),plus=await evaluate(b);
      columns.push(plus.residual.map((v,j)=>(v-minus.residual[j])/(b[i]-a[i])));
    }
    const step=leastSquaresStep(columns,best.residual).map((v,i)=>Math.max(-axes[i].maxStep,Math.min(axes[i].maxStep,v)));
    let improved=false;
    for(const scale of [1,.5,.25]) {
      const x=best.x.map((v,i)=>Math.max(axes[i].lo,Math.min(axes[i].hi,v+scale*step[i])));
      const trial=await evaluate(x);
      if(trial.score<best.score-1e-5) {best=trial;improved=true;break;}
    }
    if(!improved)break;
  }
  const after=rms(best.errors.late), front=rms(best.errors.front);
  const values=unpack(best.x);
  const atLimit=axes.some((a,i)=>Math.min(best.x[i]-a.lo,a.hi-best.x[i])<.002);
  let accepted=after<before*.9 && after<=1.5 && front<=1.5 && Math.max(...best.errors.front.map(Math.abs))<=4.5;
  // Do not accept a correction that only works for the struck random phase.
  let validation;
  if(accepted) {
    const other=(seed+911)>>>0;
    const b=await probe(baseline,other),e=await probe(edited,other),c=await probe(values,other);
    const t=holdTargets(b,e),old=t.errors(e),next=t.errors(c);
    validation={before:rms(old.late),after:rms(next.late),front:rms(next.front)};
    accepted=validation.after<=Math.max(.2,validation.before) && validation.after<=1.75 && validation.front<=1.75;
  }
  return {accepted,values:accepted?values:edited.slice(),before,after,front,atLimit,validation,evaluations,
    reason:accepted?(after>.5?'Decay partly held':'Decay held'):atLimit?'Cannot hold decay within control limits':'Correction rejected: decay or bloom changed too much'};
}
