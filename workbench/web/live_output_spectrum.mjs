import {spectrumHistogram} from './radiation_response.mjs';

// Transparent full-mix tap before browser monitor gain/limiter. Analysis does
// not delay audio. Web Audio's periodic Blackman window has gain .42 and
// ENBW (a0²+(a1²+a2²)/2)/a0² bins; correct these explicitly, never auto-level.
export class LiveOutputSpectrum {
  constructor(context) {
    this.context=context;
    this.node=new AnalyserNode(context,{fftSize:4096,smoothingTimeConstant:.4});
    this.values=new Float32Array(this.node.frequencyBinCount);
  }
  read() {
    if(this.context.state!=='running')return null;
    this.node.getFloatFrequencyData(this.values);
    const correction=20*Math.log10(2/.42);
    let peak=-Infinity;
    for(let i=0;i<this.values.length;i++) {
      this.values[i]+=correction;
      peak=Math.max(peak,this.values[i]);
    }
    if(peak < -110)return null;
    const sampleRate=this.context.sampleRate,size=this.node.fftSize;
    const noiseBandwidthHz=(.42**2+(.5**2+.08**2)/2)/.42**2*sampleRate/size;
    return spectrumHistogram({values:this.values,frames:1,bins:this.values.length,
      size,hop:size,sampleRate,noiseBandwidthHz});
  }
}
