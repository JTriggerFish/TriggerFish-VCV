// Silent, disposable CDP tab: never touch the user's playing tab or saved fits.
import assert from 'node:assert/strict';

const endpoint = 'http://127.0.0.1:9223';
const page = await (await fetch(endpoint+'/json/new?about:blank', {method:'PUT'})).json();
const socket = new WebSocket(page.webSocketDebuggerUrl);
await new Promise(resolve => socket.addEventListener('open', resolve, {once:true}));
let sequence = 0;
const pending = new Map(), errors = [];
socket.onmessage = event => {
  const m = JSON.parse(event.data);
  if (m.method === 'Runtime.exceptionThrown') errors.push(m.params.exceptionDetails.text);
  if (!pending.has(m.id)) return;
  const {resolve,reject} = pending.get(m.id); pending.delete(m.id);
  m.error ? reject(Error(m.error.message)) : resolve(m.result);
};
const call = (method,params={}) => new Promise((resolve,reject) => {
  const id = ++sequence; pending.set(id,{resolve,reject});
  socket.send(JSON.stringify({id,method,params}));
});
const evaluate = async expression => {
  const r = await call('Runtime.evaluate',{expression,awaitPromise:true,returnByValue:true});
  if (r.exceptionDetails) throw Error(JSON.stringify(r.exceptionDetails));
  return r.result.value;
};



try {
  await call('Runtime.enable');await call('Page.enable');await call('Network.enable');
  await call('Network.setCacheDisabled',{cacheDisabled:true});
  await call('Page.addScriptToEvaluateOnNewDocument',{source:
    'Object.defineProperty(navigator,"requestMIDIAccess",{value:undefined});'});
  await call('Page.navigate',{url:'http://127.0.0.1:8765/'});
  const results=await evaluate(`(async()=>{
    const check=(v,m)=>{if(!v)throw Error(m);};
    await new Promise((resolve,reject)=>{
      const deadline=performance.now()+45000;
      const poll=()=>document.querySelector('#instrument-calibration')?.options.length>2
        ?resolve(true):performance.now()>deadline?reject(Error('Not ready')):setTimeout(poll,100);poll();
    });
    const selector=document.querySelector('#instrument-calibration');
    selector.value='gong-standard';await selector.onchange();
    check(!document.querySelector('#modal-preset'),'Quick shapes remain');
    const panel=document.querySelector('#modal-templates');
    const field=k=>panel.querySelector('[data-template-key="'+k+'"]');
    const edit=(k,v)=>{field(k).value=v;field(k).dispatchEvent(new Event('input',{bubbles:true}));};
    const family=panel.querySelector('[aria-label="Modal formula"]');
    const button=panel.querySelector('.template-actions button');
    edit('count',32);family.value='membrane';family.dispatchEvent(new Event('change',{bubbles:true}));
    check(!button.disabled&&field('count').value==='32','Membrane rejects 32 modes');
    const capture=async()=>{
      let blob;const create=URL.createObjectURL,click=HTMLAnchorElement.prototype.click;
      URL.createObjectURL=b=>{blob=b;return create.call(URL,b);};
      HTMLAnchorElement.prototype.click=function(){};
      try{document.querySelector('#save-fit').click();}
      finally{URL.createObjectURL=create;HTMLAnchorElement.prototype.click=click;}
      return Object.assign({},...JSON.parse(await blob.text()).instrument.nodes.map(n=>n.parameters));
    };
    const before=await capture();button.click();const after=await capture();
    const {MembraneRatios}=await import('./modal_templates.mjs');
    for(let i=0;i<32;i++)check(after['resolved_frequency_'+i]===55*MembraneRatios[i],'Incorrect membrane mode '+i);
    for(const key of Object.keys(before))if(!key.startsWith('resolved_'))check(before[key]===after[key],'Generator changed '+key);
    check(document.querySelectorAll('#modal-editor .modal-node').length===32,'Not all modes visible');
    edit('stretch',.6);button.click();const stretched=await capture();
    for(let i=0;i<4;i++)check(stretched['resolved_frequency_'+i]===after['resolved_frequency_'+i],'Protected core moved');
    check(stretched.resolved_frequency_31>after.resolved_frequency_31,'Upper stretch ineffective');
    edit('fundamental',10000);
    check(button.disabled&&field('count').getAttribute('aria-invalid')==='true','Out-of-range count silently accepted');
    edit('fundamental',55);
    check(!button.disabled&&field('count').parentElement.querySelector('[type=range]').value==='32','Range thumb did not recover');
    edit('turbulence',0);button.click();const clean=await capture();
    for(let i=0;i<32;i++)check(clean['resolved_turbulence_'+i]===0,'Noisiness not applied');
    selector.value='kick-standard';await selector.onchange();
    const kick=document.querySelector('#kick-modal-templates [data-template-key=count]');
    check(kick.max==='16','Kick capacity changed');
    const {mountModalTemplates}=await import('./modal_template_controls.mjs');
    const scratch=document.createElement('div');
    const options={capacity:32,minimumFrequency:20,maximumFrequency:15000,apply:()=>{}};
    let oldEvents=0,newEvents=0;
    const old=mountModalTemplates(scratch,options);
    old.onPitchChange(()=>oldEvents++);old.destroy();
    const fresh=mountModalTemplates(scratch,options);
    fresh.onPitchChange(()=>newEvents++);
    scratch.dispatchEvent(new Event('change',{bubbles:true}));
    check(oldEvents===0&&newEvents===1,'Rebuilt generator retains old pitch or validation listeners');
    fresh.destroy();scratch.dispatchEvent(new Event('change',{bubbles:true}));
    check(newEvents===1,'Destroyed generator still responds');
    check(document.querySelector('#error-banner').hidden,'Unexpected error banner');
    return {membraneModes:32,exactSavedFrequencies:true,unchangedBloomAndDecay:true,stretch:true,noisiness:true,kickCapacity:16};
  })()`);
  assert.deepEqual(errors,[]);console.log(results);
} finally {
  socket.close();await fetch(endpoint+'/json/close/'+page.id);
}
