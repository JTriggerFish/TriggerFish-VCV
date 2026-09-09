// Small damped least-squares solve. All quantities are design-time coordinates;
// no inverse damping, gain matching or solver state enters the voice.
export const rms = values => Math.sqrt(values.reduce((a,x)=>a+x*x,0)/Math.max(1,values.length));

export function leastSquaresStep(columns, residual) {
  const n=columns.length;
  const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
  const matrix=columns.map(a=>columns.map(b=>dot(a,b)));
  const ridge=Math.max(1e-5,...matrix.map((r,i)=>r[i]*1e-4));
  const rhs=columns.map(c=>-dot(c,residual));
  for(let i=0;i<n;i++)matrix[i][i]+=ridge;
  for(let i=0;i<n;i++) {
    let pivot=i;
    for(let j=i+1;j<n;j++)if(Math.abs(matrix[j][i])>Math.abs(matrix[pivot][i]))pivot=j;
    [matrix[i],matrix[pivot]]=[matrix[pivot],matrix[i]];
    [rhs[i],rhs[pivot]]=[rhs[pivot],rhs[i]];
    for(let j=i+1;j<n;j++) {
      const scale=matrix[j][i]/matrix[i][i];
      for(let k=i;k<n;k++)matrix[j][k]-=scale*matrix[i][k];
      rhs[j]-=scale*rhs[i];
    }
  }
  const result=new Array(n).fill(0);
  for(let i=n-1;i>=0;i--)result[i]=(rhs[i]-matrix[i].reduce((s,x,j)=>s+(j>i?x*result[j]:0),0))/matrix[i][i];
  return result;
}

export function compensationCoordinates(descriptors, baseline, edited) {
  return descriptors.filter(d=>{
    if(d.key==='bloom_energy_acceleration')return baseline[d.index]===edited[d.index];
    const match=/^body_decay_seconds_(\d)$/.exec(d.key);
    if(!match)return false;
    const i=Number(match[1]);
    return i===0 || i===7 || edited[descriptors.find(x=>x.key===`body_decay_active_${i}`).index]>=.5;
  }).map(d=>{
    const log=d.key.startsWith('body_decay_');
    const encode=x=>log?Math.log(x):x, decode=x=>log?Math.exp(x):x;
    const origin=encode(edited[d.index]), reach=log?Math.log(2):.12;
    return {d,origin,decode,step:log?.03:.015,maxStep:log?.35:.06,
      lo:Math.max(encode(d.minimum),origin-reach),hi:Math.min(encode(d.maximum),origin+reach)};
  });
}
