// Verify target -> editable patch -> saved fit in an isolated browser tab.
const endpoint="http://127.0.0.1:9223";
const targets=process.argv.slice(2);
if(!targets.length) targets.push("snare-standard");
const page=await(await fetch(endpoint+"/json/new?about:blank",{method:"PUT"})).json();
const socket=new WebSocket(page.webSocketDebuggerUrl);
await new Promise(resolve=>socket.addEventListener("open",resolve,{once:true}));
let sequence=0;
const pending=new Map(), errors=[];
socket.onmessage=event=>{
  const m=JSON.parse(event.data);
  if(m.method==="Runtime.exceptionThrown")errors.push(m.params.exceptionDetails.text);
  if(pending.has(m.id)){
    const {resolve,reject}=pending.get(m.id);pending.delete(m.id);
    m.error?reject(Error(m.error.message)):resolve(m.result);
  }
};
const call=(method,params={})=>new Promise((resolve,reject)=>{
  const id=++sequence;pending.set(id,{resolve,reject});socket.send(JSON.stringify({id,method,params}));
});
const evaluate=async expression=>{
  const result=await call("Runtime.evaluate",{expression,awaitPromise:true,returnByValue:true});
  if(result.exceptionDetails)throw Error(JSON.stringify(result.exceptionDetails));
  return result.result.value;
};
try {
  await call("Runtime.enable"); await call("Page.enable");
  // This save/load probe must not request access to the user's MIDI hardware.
  await call("Page.addScriptToEvaluateOnNewDocument", {source: `
    Object.defineProperty(navigator,"requestMIDIAccess",{value:undefined});
    window.renderRequests=[];
    window.Worker=class extends Worker {
      constructor(url, options){super(url,options);this.url=String(url);}
      postMessage(data,...rest){
        if(this.url==='render_worker.mjs')renderRequests.push(data);
        super.postMessage(data,...rest);
      }
    };`});
  await call("Page.navigate",{url:"http://127.0.0.1:8765/"});
  await evaluate(`new Promise((resolve,reject)=>{
    const deadline=performance.now()+45000;
    const poll=()=>{
      if(document.getElementById('instrument-calibration')?.options.length>2)resolve(true);
      else if(performance.now()>deadline)reject(Error('Targets did not load'));
      else setTimeout(poll,100);
    };poll();
  })`);
  await evaluate(`(async()=>{
    const {paintLimiterMeter}=await import('./limiter_meter.mjs');
    const meter=document.getElementById('limiter-meter');
    if(!meter?.getClientRects().length)throw Error('Limiter meter is not visible');
    paintLimiterMeter(meter,{active:true,recentDb:10.7,maximumDb:10.7});
    if(!meter.classList.contains('is-limiting') ||
       meter.querySelector('meter').value!==10.7 ||
       meter.querySelector('.limiter-label').textContent!=='LIMITING')
      throw Error('Limiter gain reduction is not displayed');
    paintLimiterMeter(meter,{active:false,recentDb:0,maximumDb:10.7});
    if(!meter.classList.contains('has-limited') ||
       !meter.querySelector('button').textContent.includes('10.7'))
      throw Error('Limiter maximum was lost');
    document.getElementById('limiter-reset').click();
    if(meter.classList.contains('has-limited') ||
       !meter.querySelector('button').textContent.includes('0.0'))
      throw Error('Limiter maximum reset failed');
  })()`);
  for(const target of targets){
    const result=await evaluate(`(async()=>{
      const id=${JSON.stringify(target)};
      const {referenceCalibration}=await import('./reference_calibration_library.mjs');
      const expected=referenceCalibration(id);
      if(!expected)throw Error('No fitted calibration for '+id);
      const selector=document.getElementById('instrument-calibration');
      selector.value=id; await selector.onchange();
      if(expected.instrument.recipe==='metal.cymbal.v1') {
        for(const [key,label] of [['bloom_energy_acceleration','Concentration dependence'],
                                 ['bloom_energy_sensitivity','Energy sensitivity'],
                                 ['field_motion_depth','Amount'],
                                 ['field_motion_rate','Speed'],
                                 ['field_motion_sharing','Shimmer moves together']]) {
          const row=document.querySelector('[data-fit-key="'+key+'"]');
          if(!row?.querySelector('input[type=range]') || !row.textContent.includes(label) ||
             !row.dataset.tooltip || !row.getClientRects().length)
            throw Error('Missing visible diffusion control/help: '+key);
        }
      }
      let captured;
      const create=URL.createObjectURL,click=HTMLAnchorElement.prototype.click;
      URL.createObjectURL=blob=>{captured=blob;return create.call(URL,blob);};
      HTMLAnchorElement.prototype.click=function(){if(!this.download)click.call(this);};
      try{document.getElementById('save-fit').click();}
      finally{URL.createObjectURL=create;HTMLAnchorElement.prototype.click=click;}
      if(!captured)throw Error('No saved fit');
      const actual=JSON.parse(await captured.text());
      const values=x=>Object.assign({},...x.instrument.nodes.map(n=>n.parameters));
      const want=values(expected), got=values(actual);
      if(Object.keys(want).length!==Object.keys(got).length ||
         Object.keys(want).some(k=>want[k]!==got[k]))throw Error('Wrong parameters for '+id);
      for(const [key,value] of Object.entries(expected.controls.event))
        if(actual.controls.event[key]!==value)throw Error('Wrong event '+key+' for '+id);
      if(actual.reference.sha256!==expected.reference.sha256)throw Error('Wrong reference for '+id);
      if(actual.reference.cell?.onset_seconds!==expected.reference.cell?.onset_seconds)
        throw Error('Wrong reference onset for '+id);
      if(actual.reference.referenceGainDb!==expected.reference.referenceGainDb)
        throw Error('Wrong reference gain for '+id);
      return {target:id,parameters:Object.keys(got).length,reference:actual.reference.sha256};
    })()`);
    console.log(JSON.stringify(result));
  }
  await evaluate(`(async()=>{
    const {default:fits}=await import('./texture_trials.fit.json',{with:{type:'json'}});
    const select=document.getElementById('texture-trial');
    if(fits.length && (!select?.getClientRects().length || select.disabled))
      throw Error('Texture trials are not discoverable');
    const values=fit=>Object.assign({},...fit.instrument.nodes.map(n=>n.parameters));
    for(const fit of fits){
      select.value=fit.id; await select.onchange();
      const eq=document.querySelector('[data-fit-key="output_eq_enabled"] input');
      if(eq?.checked)throw Error('Texture trial enabled output EQ');
      const expected=values(fit);
      if(expected.output_eq_enabled!==0)throw Error('Trial is not EQ-free');
      let captured;
      const create=URL.createObjectURL,click=HTMLAnchorElement.prototype.click;
      URL.createObjectURL=blob=>{captured=blob;return create.call(URL,blob);};
      HTMLAnchorElement.prototype.click=function(){if(!this.download)click.call(this);};
      try{document.getElementById('save-fit').click();}
      finally{URL.createObjectURL=create;HTMLAnchorElement.prototype.click=click;}
      const actual=JSON.parse(await captured.text()), got=values(actual);
      if(Object.keys(expected).some(key=>expected[key]!==got[key]))
        throw Error('Texture trial did not restore its controls');
      if(actual.reference.sha256!==fit.reference.sha256 ||
         actual.reference.referenceGainDb!==fit.reference.referenceGainDb)
        throw Error('Texture trial changed reference or gain');
      for(const [key,value] of Object.entries(fit.controls.event))
        if(actual.controls.event[key]!==value)
          throw Error('Texture trial changed saved gesture: '+key);
      if(document.getElementById('snapshot-name').value!==fit.name)
        throw Error('Texture trial name was not restored');
    }
    return true;
  })()`);
  await evaluate(`(async()=>{
    const wait=async test=>{const deadline=performance.now()+30000;while(!test()){
      if(performance.now()>deadline)throw Error('Snapshot test timed out');
      await new Promise(r=>setTimeout(r,20));}};
    const ready=()=>document.getElementById('status').textContent==='Ready';
    await wait(ready);
    const snapshot=name=>{document.getElementById('snapshot-name').value=name;
      document.getElementById('snapshot').click();
      return [...document.querySelectorAll('.snapshot-chip')].at(-1);};
    const cached=snapshot('Completed audio');
    const input=document.querySelector('[data-fit-key=body_brightness] input[type=range]');
    const changed=.5; // normalized slider position: -24 dB/oct
    input.value=changed;input.dispatchEvent(new Event('input',{bubbles:true}));
    const pending=snapshot('Pending audio');
    await wait(ready);
    // Restoring a pending snapshot must render its controls, not reuse old PCM.
    cached.click();await wait(ready);renderRequests.length=0;
    pending.click();await wait(()=>renderRequests.some(r=>!r.cancel));await wait(ready);
    // A cached restore must cancel a newer render which could otherwise replace it.
    const newer=document.querySelector('[data-fit-key=body_brightness] input[type=range]');
    newer.value=.55;newer.dispatchEvent(new Event('input',{bubbles:true}));
    await wait(()=>document.getElementById('render-time').textContent.startsWith('Rendering'));
    renderRequests.length=0;cached.click();await wait(ready);
    if(!renderRequests.some(r=>r.cancel))throw Error('Cached restore did not cancel old render');
    // Save fit must reflect visible edits, not the selected snapshot's old values.
    const visible=document.querySelector('[data-fit-key=body_brightness] input[type=range]');
    visible.value=changed;visible.dispatchEvent(new Event('input',{bubbles:true}));
    let blob;const create=URL.createObjectURL,click=HTMLAnchorElement.prototype.click;
    URL.createObjectURL=value=>{blob=value;return create.call(URL,value);};
    HTMLAnchorElement.prototype.click=function(){if(!this.download)click.call(this);};
    try{document.getElementById('save-fit').click();}
    finally{URL.createObjectURL=create;HTMLAnchorElement.prototype.click=click;}
    const fit=JSON.parse(await blob.text());
    if(fit.instrument.nodes.find(n=>n.id==='body').parameters.body_brightness!==-24)
      throw Error('Save fit lost edits after snapshot selection');
  })()`);
  if(errors.length)throw Error(errors.join('\n'));
} finally {
  socket.close();await fetch(endpoint+"/json/close/"+page.id);
}
