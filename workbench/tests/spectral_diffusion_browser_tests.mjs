// Optional silent integration check in a disposable CDP tab, never the user's.
import {writeFile} from "node:fs/promises";
const endpoint = process.env.TF_CDP_URL ?? "http://127.0.0.1:9223";
const page = await (await fetch(`${endpoint}/json/new?about:blank`, {method:"PUT"})).json();
const socket = new WebSocket(page.webSocketDebuggerUrl);
await new Promise(resolve => socket.addEventListener("open", resolve, {once:true}));
let sequence = 0;
const pending = new Map();
socket.onmessage = event => {
  const message = JSON.parse(event.data), item = pending.get(message.id);
  if (!item) return;
  pending.delete(message.id);
  message.error ? item.reject(Error(message.error.message)) : item.resolve(message.result);
};
const call = (method, params={}) => new Promise((resolve,reject) => {
  const id=++sequence; pending.set(id,{resolve,reject});
  socket.send(JSON.stringify({id,method,params}));
});
try {
  await call("Page.enable");
  await call("Network.enable");
  await call("Network.setCacheDisabled", {cacheDisabled:true});
  await call("Emulation.setDeviceMetricsOverride",{width:1800,height:1200,deviceScaleFactor:1,mobile:false});
  // This silent generator test must not request access to the user's MIDI devices.
  await call("Page.addScriptToEvaluateOnNewDocument", {source:
    `Object.defineProperty(navigator,'requestMIDIAccess',{value:undefined,configurable:true});
     window.testWorkers=[];
     window.Worker=class extends Worker { constructor(url,options) {
       super(url,options);window.testWorkers.push({url:String(url),worker:this});
     }};`});
  await call("Page.navigate", {url:"http://127.0.0.1:8765"});
  await new Promise(resolve => setTimeout(resolve, 1000));
  const result = await call("Runtime.evaluate", {awaitPromise:true,returnByValue:true,
    expression:`(async () => {
      const wait = async predicate => {
        for(let i=0;i<300;++i) {
          if(predicate()) return;
          await new Promise(resolve=>setTimeout(resolve,100));
        }
        throw Error('Timed out: '+document.getElementById('status')?.textContent);
      };
      const row = key => document.querySelector('[data-fit-key="'+key+'"]');
      const ready = () => document.getElementById('status')?.textContent.startsWith('Ready');
      await wait(ready);
      const results=[];
      for (const id of ['crash-standard','gong-standard','ride-standard','hihat-standard']) {
        const select=document.getElementById('instrument-calibration');
        if(![...select.options].some(o=>o.value===id)) throw Error('Missing '+id);
        select.value=id;select.dispatchEvent(new Event('change',{bubbles:true}));
        await new Promise(resolve=>setTimeout(resolve,200));await wait(ready);
        if(!document.getElementById('error-banner').hidden)
          throw Error('Preset load reported: '+document.getElementById('error-message').textContent);
        for (const obsolete of ['bloom_spectral_diffusion','field_relaxed_turbulence',
            'bloom_phase_diffusion','field_exchange'])
          if(row(obsolete)) throw Error('Obsolete control visible: '+obsolete);
        if(!row('bloom_rate').textContent.includes('Diffusion strength') ||
           row('bloom_rate').textContent.includes('oct/s') ||
           !row('field_turbulence').textContent.includes('Packet noisiness'))
          throw Error('Wrong model labels');
        results.push(id);
      }
      // Cancelled preset fetches must not overwrite the current gesture or
      // produce a false unavailable-reference banner.
      const presets=document.getElementById('instrument-calibration');
      for(const id of ['gong-standard','ride-standard','crash-standard']) {
        presets.value=id;presets.dispatchEvent(new Event('change'));
      }
      await new Promise(resolve=>setTimeout(resolve,200));await wait(ready);
      if(!document.getElementById('error-banner').hidden)
        throw Error('Cancelled preset load reported an error');
      // Drop one response, then inject a decode failure while a job is pending.
      // A subsequent edit must render and analyse again, with no audio device.
      for(const name of ['render_worker.mjs','analysis_worker.mjs']) {
        const target=window.testWorkers.find(item=>item.url===name).worker;
        const post=target.postMessage.bind(target);let dropped=false;
        target.postMessage=()=>{dropped=true;};
        presets.value='ride-standard';presets.dispatchEvent(new Event('change'));
        await wait(()=>dropped);
        target.dispatchEvent(new MessageEvent('messageerror'));
        if(document.getElementById('error-banner').hidden)
          throw Error('Worker decode failure was swallowed: '+name+' handler='+
            String(target.onmessageerror)+' status='+document.getElementById('status').textContent);
        target.postMessage=post;
        document.getElementById('dismiss-error').click();
        presets.value='crash-standard';presets.dispatchEvent(new Event('change'));
        await new Promise(resolve=>setTimeout(resolve,200));await wait(ready);
      }
      const {PercussionEngine}=await import('./engine.mjs');
      const {FitControls}=await import('./fit_controls.mjs');
      const engine=await PercussionEngine.create(44100,0);
      const controls=new FitControls({descriptors:engine.parameters,
        state:{macros:engine.parameters.map(d=>d.defaultValue)},onChange:()=>{}});
      document.getElementById('modal-editor').replaceChildren();
      controls.buildResolvedEditor();
      const host=document.getElementById('modal-templates');
      if(!host.querySelector('details').open) throw Error('Generator is hidden');
      if(!host.querySelector('details').open ||
         host.querySelector('[aria-label="Modal formula"]').value!=='harmonic')
        throw Error('Harmonic generator not visible by default');
      const note=host.querySelector('[aria-label="Generator base note"]');
      const octave=host.querySelector('[aria-label="Generator base octave"]');
      note.value='11';note.dispatchEvent(new Event('change'));
      octave.value='2';octave.dispatchEvent(new Event('change'));
      const hz=host.querySelector('[data-template-key="fundamental"]');
      if(Math.abs(Number(hz.value)-123.470825314)>.000001) throw Error('Wrong note conversion');
      host.querySelector('[data-template-key="count"]').value=3;
      host.querySelector('button').click();
      const parameters=Object.fromEntries(engine.parameters.map(d=>[d.key,controls.state.macros[d.index]]));
      for(let i=0;i<3;++i)
        if(Math.abs(parameters['resolved_frequency_'+i]-Number(hz.value)*(i+1))>.00001 ||
           parameters['resolved_turbulence_'+i]!==1) throw Error('Incorrect harmonic modes/noisiness disabled');
      if(parameters.resolved_level_3!==-72) throw Error('Stale ungenerated mode');
      hz.value='127.3';hz.dispatchEvent(new Event('input'));
      if(note.value!=='custom' || !octave.disabled) throw Error('Custom frequency not preserved');
      hz.dispatchEvent(new MouseEvent('dblclick'));
      if(hz.value!=='55' || note.value!=='9' || octave.value!=='1') throw Error('Reset desynchronized note');
      const count=host.querySelector('[data-template-key="count"]');
      const generate=host.querySelector('button');
      count.value=32;count.dispatchEvent(new Event('input',{bubbles:true}));generate.click();
      if(document.querySelectorAll('#modal-editor .modal-bar').length!==32)
        throw Error('32 modes failed: '+host.querySelector('.template-status').textContent+
          ' count='+count.value+' bars='+document.querySelectorAll('#modal-editor .modal-bar').length);
      const before=JSON.stringify(controls.state.macros);
      count.value=33;count.dispatchEvent(new Event('change',{bubbles:true}));
      if(!generate.disabled || !host.querySelector('.template-status').textContent.includes('1–32') ||
         document.getElementById('error-banner').hidden || JSON.stringify(controls.state.macros)!==before)
        throw Error('Invalid count not clearly rejected without changing patch');
      count.value=16;count.dispatchEvent(new Event('input',{bubbles:true}));
      const formula=host.querySelector('[aria-label="Modal formula"]');
      formula.value='membrane';formula.dispatchEvent(new Event('change',{bubbles:true}));
      if(count.max!=='16' || generate.disabled) throw Error('Membrane limit incorrect');
      count.value=17;count.dispatchEvent(new Event('input',{bubbles:true}));
      if(!generate.disabled) throw Error('Membrane silently truncated');
      formula.value='harmonic';formula.dispatchEvent(new Event('change'));
      hz.value=1000;hz.dispatchEvent(new Event('input',{bubbles:true}));
      if(count.max!=='15' || !generate.disabled) throw Error('Frequency ceiling not enforced');
      hz.value=110;hz.dispatchEvent(new Event('input',{bubbles:true}));
      count.value=8;count.dispatchEvent(new Event('input',{bubbles:true}));generate.click();
      const values=controls.state.macros.slice();
      const set=(key,value)=>{values[engine.parameters.find(d=>d.key===key).index]=value;};
      set('bloom_rate',0);set('direct_gain',0);set('field_turbulence_slope',0);
      set('field_phase_bandwidth',1);set('field_packet_spread',1);
      const render=noise=>{set('field_turbulence',noise);engine.reset();
        return engine.render({seconds:.25,parameters:values,seed:47});};
      const pure=render(0), noisy=render(1);
      const energy=pure.reduce((sum,x)=>sum+x*x,0);
      const delta=pure.reduce((sum,x,i)=>sum+(x-noisy[i])**2,0);
      if(!noisy.every(Number.isFinite) || energy<=0 || delta/energy<.01)
        throw Error('Generated modes still do not respond to noisiness in DSP');
      const guide=document.getElementById('harmonic-guide');
      guide.checked=true;guide.dispatchEvent(new Event('change'));
      if(controls.resolvedEditor.harmonicGuide.fundamentalHz!==110) throw Error('Guide has a different base pitch');
      hz.value=10;hz.dispatchEvent(new Event('input',{bubbles:true}));
      count.value=4;count.dispatchEvent(new Event('input',{bubbles:true}));generate.click();
      if(generate.disabled || controls.resolvedEditor.point(0).frequency!==10 ||
         controls.resolvedEditor.harmonicGuide.fundamentalHz!==10)
        throw Error('Sub-bass generation/editor/guide failed');
      if(controls.resolvedEditor.options.minimumFrequency!==20 ||
         controls.resolvedEditor.svg.querySelectorAll('.modal-bar').length!==3)
        throw Error('Display must start at 20 Hz without stacking lower modes at the edge');
      engine.setConfiguration(controls.state.macros);engine.reset();engine.trigger({seed:57});
      const lowAudio=new Float32Array(44100);engine.processTo(lowAudio,0,lowAudio.length);
      if(!lowAudio.every(Number.isFinite) || !lowAudio.some(x=>Math.abs(x)>1e-7))
        throw Error('Sub-bass patch did not render through actual API');
      const stretch=host.querySelector('[data-template-key="stretch"]');
      const core=host.querySelector('[data-template-key="harmonicCore"]');
      if(!core || Number(core.value)!==4) throw Error('Missing default harmonic core');
      count.value=8;count.dispatchEvent(new Event('input',{bubbles:true}));
      stretch.value=.5;stretch.dispatchEvent(new Event('input',{bubbles:true}));generate.click();
      if(controls.resolvedEditor.point(0).frequency!==10 ||
         controls.resolvedEditor.point(3).frequency!==40 ||
         Math.abs(controls.resolvedEditor.point(7).frequency-80*Math.sqrt(1.25))>.0001)
        throw Error('Stretch failed to protect the low core or expand higher modes');
      core.value=2;core.dispatchEvent(new Event('input',{bubbles:true}));generate.click();
      if(controls.resolvedEditor.point(1).frequency!==20 ||
         !(controls.resolvedEditor.point(3).frequency>40))
        throw Error('Harmonic core control did not change the stretch onset');
      const {reportError}=await import('./error_banner.mjs');
      reportError(new Error('Persistent test error'));
      document.getElementById('status').textContent='Ready';
      if(document.getElementById('error-banner').hidden) throw Error('Ready cleared error');
      window.dispatchEvent(new ErrorEvent('error',{message:'Unhandled test error'}));
      if(!document.getElementById('error-message').textContent.includes('Unhandled test')) throw Error('Uncaught error swallowed');
      window.dispatchEvent(new PromiseRejectionEvent('unhandledrejection',{reason:new Error('Async test error'),promise:Promise.resolve()}));
      if(!document.getElementById('error-message').textContent.includes('Async test')) throw Error('Async error swallowed');
      document.getElementById('dismiss-error').click();
      if(!document.getElementById('error-banner').hidden) throw Error('Dismiss failed');
      reportError(new Error('Example validation notice for visual review — this disposable tab only.'), 'Mode generator');
      host.scrollIntoView({block:'center'});
      controls.resolvedEditor.resizeObserver.disconnect();engine.destroy();
      return {ready:results,newModelOnly:true,harmonicGenerator:true,limits:true,
        noiseDifference:delta/energy,errorBanner:true,audioPlayed:false};
    })()`});
  if(result.exceptionDetails) throw Error(JSON.stringify(result.exceptionDetails));
  console.log(JSON.stringify(result.result.value));
  const shot=await call("Page.captureScreenshot",{format:"png"});
  await writeFile("build/modal-generator-review.png",Buffer.from(shot.data,"base64"));
} catch (error) {
  console.error(error); process.exitCode=1;
} finally {
  const closed=new Promise(resolve=>socket.addEventListener("close",resolve,{once:true}));
  socket.close();await closed;await fetch(endpoint+"/json/close/"+page.id);
}
