import {fft, windowSamples} from "./analysis.mjs";

// Fixed analysis coordinates, independent of the display FFT and colour scale.
const Regions = [[0,.05],[.05,.1],[.1,.2],[.2,.3],[.3,.45],
  [1,1.5],[1.5,2],[2,3],[3,4],[4,5],[5,6]];
export const FrontCells = 24 * 5;

export function decayEnvelope(samples, sampleRate) {
  if (!Number.isFinite(sampleRate) || sampleRate < 8000 ||
      samples.length < Math.round(6 * sampleRate)) throw Error("Hold decay needs a six-second render");
  const size=4096, hop=Math.round(sampleRate*.03), window=windowSamples("hann",size);
  const real=new Float64Array(size), imaginary=new Float64Array(size);
  const bands=Array.from({length:24},(_,i)=>[80*(Math.min(16000,.45*sampleRate)/80)**(i/24),
    80*(Math.min(16000,.45*sampleRate)/80)**((i+1)/24)]);
  const binBand=Array.from({length:size/2+1},(_,i)=>bands.findIndex(([lo,hi])=>i*sampleRate/size>=lo && i*sampleRate/size<hi));
  const power=Regions.map(()=>new Float64Array(24)), counts=new Uint32Array(Regions.length);
  for(let centre=0;centre<samples.length;centre+=hop) {
    const region=Regions.findIndex(([lo,hi])=>centre/sampleRate>=lo && centre/sampleRate<hi);
    if(region<0)continue;
    for(let i=0;i<size;i++) {
      const value=samples[centre+i-size/2] ?? 0;
      if(!Number.isFinite(value))throw Error("Nonfinite hold-decay render");
      real[i]=value*window[i];imaginary[i]=0;
    }
    fft(real,imaginary);counts[region]++;
    for(let bin=1;bin<binBand.length;bin++)if(binBand[bin]>=0)
      power[region][binBand[bin]]+=(real[bin]**2+imaginary[bin]**2)/(size*size);
  }
  return Float64Array.from(power.flatMap((row,i)=>Array.from(row,x=>x/Math.max(1,counts[i]))));
}

export function holdTargets(baseline, edited) {
  if(baseline.length!==264 || edited.length!==264 ||
      [...baseline,...edited].some(x=>!Number.isFinite(x)||x<0))throw Error("Invalid decay measurements");
  const peak=Math.max(...baseline);
  if(peak<1e-18)throw Error("Hold decay needs an audible starting sound");
  const floor=peak*1e-7;
  const db=x=>10*Math.log10(Math.max(floor,x));
  const late=[], front=[];
  baseline.forEach((value,index)=>{
    if(index>=FrontCells && value>peak*10**(-55/10))late.push(index);
    if(index<FrontCells && Math.max(value,edited[index])>peak*1e-5)front.push(index);
  });
  if(late.length<8)throw Error("Not enough audible tail to hold");
  return {
    errors: measured=> {
      if (measured.length!==264 || measured.some(x=>!Number.isFinite(x)||x<0))
        throw Error('Invalid rendered decay measurements');
      return {late:late.map(i=>db(measured[i])-db(baseline[i])),
        front:front.map(i=>db(measured[i])-db(edited[i]))};
    },
    cells:late.length,
  };
}
