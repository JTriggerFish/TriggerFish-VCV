import {radiationCurves,spectrumHistogram,eqResponseSampleRate} from './radiation_response.mjs';
const ns='http://www.w3.org/2000/svg', width=300,height=164;
const clamp=(v,a,b)=>Math.max(a,Math.min(b,v));
const x=f=>26+262*Math.log(clamp(f,20,20000)/20)/Math.log(1000);
const y=db=>16+(18-clamp(db,-36,18))*124/54;
const path=(frequencies,db)=>frequencies.map((f,i)=>`${i?'L':'M'}${x(f)},${y(db[i])}`).join(' ');
const colours=['#dcb66c','#c78ac0','#74b8d6'];
function element(tag,attributes,text) {
  const e=document.createElementNS(ns,tag);
  for(const [key,v]of Object.entries(attributes))e.setAttribute(key,v);
  if(text)e.textContent=text;
  return e;
}

// Graph and editable numeric fields share the existing four DSP parameters.
export class MiniEqEditor {
  constructor(host,{prefix,read,set,descriptor,state}) {
    Object.assign(this,{prefix,read,set,descriptor,state});
    this.keys=['low_cut','colour_frequency','colour_gain','high_cut'].map(k=>`${prefix}_${k}`);
    this.svg=element('svg',{viewBox:`0 0 ${width} ${height}`,class:'mini-eq',
      role:'group','aria-label':'Final output EQ: drag nodes or edit values below'});
    host.append(this.svg);
    this.fields=document.createElement('div');this.fields.className='mini-eq-values';host.append(this.fields);
    this.inputs=this.keys.map((key,i)=>this.field(key,['High-pass Hz','Colour Hz','Colour dB','Low-pass Hz'][i]));
    this.legend=document.createElement('p');this.legend.className='control-help';
    host.append(this.legend);
    this.svg.onpointerdown=e=>this.start(e);
    this.svg.onpointermove=e=>this.move(e);
    this.svg.onpointerup=this.svg.onpointercancel=this.svg.onlostpointercapture=()=>{this.drag=null;};
    this.svg.ondblclick=e=>this.reset(e);
    this.refresh();
  }

  field(key,label) {
    const wrapper=document.createElement('label');wrapper.dataset.fitKey=key;
    wrapper.textContent=label;
    const input=document.createElement('input'),d=this.descriptor(key);
    Object.assign(input,{type:'number',min:d.minimum,max:d.maximum,step:key.endsWith('gain')?.1:1});
    input.onchange=()=>{if(Number.isFinite(input.valueAsNumber))this.set(key,clamp(input.valueAsNumber,d.minimum,d.maximum));this.refresh();};
    input.ondblclick=()=>{this.set(key,d.defaultValue);this.refresh();};
    wrapper.append(input);this.fields.append(wrapper);return input;
  }

  parameters() {
    const [low,frequency,gain,high]=this.keys.map(this.read);
    return {low,frequency,gain,high,enabled:this.read(`${this.prefix}_eq_enabled`)>=.5};
  }

  refresh() {
    this.inputs.forEach((input,i)=>{if(document.activeElement!==input)
      input.value=Number(this.read(this.keys[i]).toFixed(i===2?1:1));});
    this.svg.replaceChildren();
    this.backgroundLayer=element('g',{});
    this.svg.append(this.backgroundLayer);
    this.bars=null;
    const p=this.parameters(), rate=eqResponseSampleRate(this.state);
    this.responseSampleRate=rate;
    this.background();
    for(const db of [-24,-12,0,12]) {
      this.svg.append(element('path',{d:`M26,${y(db)}H288`,class:'eq-grid'}));
      this.svg.append(element('text',{x:2,y:y(db)+3,class:'eq-label'},String(db)));
    }
    for(const f of [100,1000,10000]) {
      this.svg.append(element('path',{d:`M${x(f)},16V140`,class:'eq-grid'}));
      this.svg.append(element('text',{x:x(f),y:156,'text-anchor':'middle',class:'eq-label'},f===100?'100':`${f/1000}k`));
    }
    const maximum=Math.min(20000,.499*rate);
    const frequencies=Array.from({length:160},(_,i)=>20*(maximum/20)**(i/159));
    const response=radiationCurves(p,rate,frequencies);
    response.parts.forEach((db,i)=>this.svg.append(element('path',{
      d:path(frequencies,db),stroke:colours[i],class:'eq-part',opacity:p.enabled?.6:.2})));
    this.svg.append(element('path',{d:path(frequencies,response.total),class:'eq-total'}));
    [[p.low,0],[p.frequency,p.gain],[p.high,0]].forEach(([f,db],i)=>{
      const node=element('circle',{cx:x(f),cy:y(db),r:6,fill:colours[i],
        'data-node':i,opacity:p.enabled?1:.35});
      node.append(element('title',{},['High-pass: drag left/right','Colour: drag pitch and gain','Low-pass: drag left/right'][i]));
      this.svg.append(node);
    });
  }

  background() {
    if(this.responseSampleRate!==eqResponseSampleRate(this.state)) {
      this.refresh(); return;
    }
    this.legend.textContent=this.state.liveEqHistogram
      ? 'Gold: reference · blue: LIVE full mix, before master/limiter. Curves: final EQ.'
      : 'Gold: reference · blue: rendered full synth · first 1 s. Curves: final EQ.';
    const spectra=[this.state.referenceSpectrum,this.state.synthesisSpectrum];
    this.histograms??=new WeakMap();
    const rows=spectra.map(s=>{
      if(!s)return null;
      if(!this.histograms.has(s))this.histograms.set(s,spectrumHistogram(s));
      return this.histograms.get(s);
    });
    if(this.state.liveEqHistogram)rows[1]=this.state.liveEqHistogram;
    const ceiling=rows[0]?Math.max(...rows[0]):-20;
    if(!this.bars)this.bars=[0,1].map(j=>Array.from({length:64},(_,i)=>{
      const rect=element('rect',{x:26+i*262/64,width:262/64,
        fill:j?'#74b8d6':'#dcb66c',opacity:.14});
      this.backgroundLayer.append(rect);return rect;
    }));
    this.bars.forEach((bars,j)=>bars.forEach((rect,i)=>{
      const db=rows[j]?.[i]??-160;
      const bar=124*clamp((db-ceiling+60)/60,0,1);
      rect.setAttribute('y',140-bar);rect.setAttribute('height',bar);
    }));
  }

  start(event) {
    const node=event.target.closest('[data-node]');
    if(!node||!this.parameters().enabled||event.button!==0)return;
    event.preventDefault();this.drag=Number(node.dataset.node);
    this.svg.setPointerCapture(event.pointerId);
  }

  move(event) {
    if(this.drag==null)return;
    const point=new DOMPoint(event.clientX,event.clientY).matrixTransform(this.svg.getScreenCTM().inverse());
    const index=[0,1,3][this.drag],key=this.keys[index],d=this.descriptor(key);
    this.set(key,clamp(20*1000**((point.x-26)/262),d.minimum,d.maximum));
    if(this.drag===1)this.set(this.keys[2],clamp(18-(point.y-16)*54/124,-18,18));
    this.refresh();
  }

  reset(event) {
    const node=event.target.closest('[data-node]');if(!node)return;
    const i=Number(node.dataset.node),keys=i===1?[this.keys[1],this.keys[2]]:[this.keys[i===0?0:3]];
    keys.forEach(key=>this.set(key,this.descriptor(key).defaultValue));this.refresh();
  }
}
