// Display-only RBJ response, matching percussion/biquad_design.hpp and the
// fixed Butterworth cut filters / Q=.8 colour peak in crash_cymbal_parameters.
const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
export const eqResponseSampleRate = state =>
  (state.liveEqHistogram && state.liveEqSampleRate) ||
  state.synthesisSpectrum?.sampleRate || state.reference?.sampleRate || 48000;
export function radiationCoefficients(type, frequency, gain, rate) {
  const w=2*Math.PI*clamp(frequency,1,.49*rate)/rate;
  const c=Math.cos(w), a=Math.sin(w)/(2*(type==='peak'?.8:Math.SQRT1_2));
  if(type==='peak') {
    const A=10**(clamp(gain,-36,36)/40);
    return [1+a*A,-2*c,1-a*A,1+a/A,-2*c,1-a/A];
  }
  const sign=type==='highpass'?1:-1, b=.5*(1+sign*c);
  return [b,-2*sign*b,b,1+a,-2*c,1-a];
}
export function responseDb(coefficients, frequency, rate) {
  const w=2*Math.PI*frequency/rate;
  const power=(a,b,c)=>(a+b*Math.cos(w)+c*Math.cos(2*w))**2+
    (b*Math.sin(w)+c*Math.sin(2*w))**2;
  return 10*Math.log10(Math.max(1e-30,power(...coefficients.slice(0,3)))/
    Math.max(1e-30,power(...coefficients.slice(3))));
}
export function radiationCurves(p, rate, frequencies) {
  const stages=[radiationCoefficients('highpass',p.low,0,rate),
    radiationCoefficients('peak',p.frequency,p.gain,rate),
    radiationCoefficients('lowpass',p.high,0,rate)];
  const parts=stages.map(c=>frequencies.map(f=>responseDb(c,f,rate)));
  return {parts,total:frequencies.map((_,i)=>p.enabled?parts.reduce((s,a)=>s+a[i],0):0)};
}

// A full-output/reference preview, NOT an invented pre-EQ contact/body signal.
// Average power over the first second; max 128 frames for bounded UI work.
export function spectrumHistogram(spectrum, count=64) {
  if(!spectrum)return null;
  const {values,frames,bins,size,hop,sampleRate}=spectrum;
  const end=Math.min(frames,Math.ceil(sampleRate/hop));
  const stride=Math.max(1,Math.ceil(end/128)), power=new Float64Array(count), n=new Uint32Array(count);
  for(let bin=1;bin<bins;bin++) {
    const f=bin*sampleRate/size;
    if(f<20||f>20000)continue;
    const bucket=Math.min(count-1,Math.floor(Math.log(f/20)/Math.log(1000)*count));
    for(let frame=0;frame<end;frame+=stride) {
      power[bucket]+=10**(values[frame*bins+bin]/10);n[bucket]++;
    }
  }
  // One-sided power density makes broadband levels independent of FFT length
  // and window. The STFT carries its measured equivalent noise bandwidth.
  const bandwidth=spectrum.noiseBandwidthHz??1.5*sampleRate/size;
  return Array.from(power,(p,i)=>10*Math.log10(Math.max(1e-16,p/Math.max(n[i],1)/(2*bandwidth))));
}
