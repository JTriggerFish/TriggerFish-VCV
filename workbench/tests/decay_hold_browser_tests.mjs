// Optional integration check: disposable silent tab, no user state or downloads.
const endpoint=process.env.TF_CDP_URL??'http://127.0.0.1:9223';
const page=await(await fetch(`${endpoint}/json/new?about:blank`,{method:'PUT'})).json();
const socket=new WebSocket(page.webSocketDebuggerUrl);
await new Promise(resolve=>socket.addEventListener('open',resolve,{once:true}));
let sequence=0;
const pending=new Map();
socket.onmessage=event=>{
  const message=JSON.parse(event.data), item=pending.get(message.id);
  if(!item)return;pending.delete(message.id);
  message.error?item.reject(Error(message.error.message)):item.resolve(message.result);
};
const call=(method,params={})=>new Promise((resolve,reject)=>{
  const id=++sequence;pending.set(id,{resolve,reject});socket.send(JSON.stringify({id,method,params}));
});
try {
  await call('Page.enable');
  await call('Network.enable');
  await call('Network.setCacheDisabled',{cacheDisabled:true});
  await call('Page.addScriptToEvaluateOnNewDocument',{source:`
    Object.defineProperty(navigator,'requestMIDIAccess',{value:undefined,configurable:true});
    window.holdResults=[]; window.holdStarted=0;
    window.Worker=class extends Worker { constructor(url,options) {
      super(url,options);
      if(String(url).includes('decay_hold_worker')) {
        window.holdStarted++;
        this.addEventListener('message',e=>{if(e.data.result)window.holdResults.push(e.data.result);});
      }
    }};
    const create=URL.createObjectURL.bind(URL);
    URL.createObjectURL=blob=>{if(blob.type==='application/json')window.savedFit=blob;return create(blob);};
    HTMLAnchorElement.prototype.click=function(){};
  `});
  await call('Page.navigate',{url:'http://127.0.0.1:8765'});
  await new Promise(resolve=>setTimeout(resolve,1000));
  const output=await call('Runtime.evaluate',{awaitPromise:true,returnByValue:true,expression:`(async()=>{
    const wait=async predicate=>{for(let i=0;i<900;i++){if(predicate())return;await new Promise(r=>setTimeout(r,100));}throw Error('Timed out');};
    const check=(value,message)=>{if(!value)throw Error(message);};
    const byId=id=>document.getElementById(id);
    const ready=()=>byId('status')?.textContent.startsWith('Ready');
    await wait(ready);
    const select=byId('instrument-calibration');
    const calibration=(await import('./reference_calibration_library.mjs')).referenceCalibration('crash-standard');
    check([...select.options].some(o=>o.value==='crash-standard'&&o.textContent===calibration.name),'Refined crash not discoverable');
    select.value='crash-standard'; select.dispatchEvent(new Event('change',{bubbles:true}));
    await new Promise(r=>setTimeout(r,500));await wait(ready);
    check(byId('snapshot-name').value===calibration.name,'Fit name not exposed');
    const row=k=>document.querySelector('[data-fit-key="'+k+'"]');
    const rate=row('bloom_rate').querySelector('input');
    const edit=(input,amount)=>{
      input.dispatchEvent(new PointerEvent('pointerdown',{bubbles:true}));
      input.value=Math.min(1,Number(input.value)+amount);
      input.dispatchEvent(new Event('input',{bubbles:true}));
      input.dispatchEvent(new Event('change',{bubbles:true}));
    };
    const beforeLevel=row('model_level_db').querySelector('output').textContent;
    edit(rate,.02);
    const requested=row('bloom_rate').querySelector('output').textContent;
    await wait(()=>window.holdResults.length===1);
    const result=window.holdResults[0];
    check(result.accepted,'Actual correction rejected: '+JSON.stringify(result));
    check(row('bloom_rate').querySelector('output').textContent===requested,'Requested bloom changed');
    check(row('model_level_db').querySelector('output').textContent===beforeLevel,'Model gain changed');
    check(document.querySelector('.decay-hold [role=status]').textContent.includes('visible controls updated'),
      'Missing result feedback: '+document.querySelector('.decay-hold [role=status]').textContent+' / '+byId('error-message').textContent);
    byId('save-fit').click();await wait(()=>window.savedFit);
    const saved=JSON.parse(await window.savedFit.text());
    const p=Object.assign({},...saved.instrument.nodes.map(n=>n.parameters));
    const {PercussionEngine}=await import('./engine.mjs');
    const engine=await PercussionEngine.create(saved.reference.sampleRate,0);
    for(const d of engine.parameters)check(p[d.key]===result.values[d.index],'Saved values differ: '+d.key);
    engine.destroy();
    edit(rate,.015);await wait(()=>window.holdStarted===2);
    document.querySelector('.decay-hold button').click();
    check(document.querySelector('.decay-hold [role=status]').textContent.includes('Cancelled'),'Cancel feedback missing');
    const count=window.holdStarted;
    edit(rate,.01);rate.dispatchEvent(new MouseEvent('dblclick',{bubbles:true}));
    await new Promise(r=>setTimeout(r,500));check(window.holdStarted===count,'Reset started correction');
    edit(rate,.01);await wait(()=>window.holdStarted===count+1);
    const volume=row('model_level_db').querySelector('input');
    volume.value=Number(volume.value)-.02;volume.dispatchEvent(new Event('input',{bubbles:true}));
    await wait(()=>document.querySelector('.decay-hold [role=status]').textContent.includes('Cancelled'));
    check(window.holdResults.length===1,'Stale job completed after edit');
    check(byId('error-banner').hidden,'Unexpected UI error');
    return {result,saveExact:true,cancel:true,reset:true,stale:true,crashDiscoverable:true};
  })()`});
  if(output.exceptionDetails)throw Error(output.exceptionDetails.exception?.description??JSON.stringify(output.exceptionDetails));
  console.log(JSON.stringify(output.result.value));
} finally {
  socket.close();await fetch(`${endpoint}/json/close/${page.id}`);
}
