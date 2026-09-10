// Silent, disposable CDP tab: never touch the user's playing tab or saved fits.
import assert from 'node:assert/strict';
import {mkdir, writeFile} from 'node:fs/promises';

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
  await call('Runtime.enable'); await call('Page.enable');
  await call('Network.enable'); await call('Network.setCacheDisabled',{cacheDisabled:true});
  await call('Emulation.setDeviceMetricsOverride',{width:1500,height:1100,deviceScaleFactor:1,mobile:false});
  await call('Page.addScriptToEvaluateOnNewDocument',{source:
    'Object.defineProperty(navigator,"requestMIDIAccess",{value:undefined});'});
  await call('Page.navigate',{url:'http://127.0.0.1:8765/'});
  await evaluate(`new Promise((resolve,reject)=>{
    const deadline=performance.now()+45000;
    const poll=()=>document.querySelector('#instrument-calibration')?.options.length>2
      ?resolve(true):performance.now()>deadline?reject(Error('Workbench did not load')):setTimeout(poll,100);
    poll();})`);
  const result = await evaluate(`(async()=>{
    const check=(v,message)=>{if(!v)throw Error(message);};
    const selector=document.querySelector('#instrument-calibration');
    selector.value='gong-standard'; await selector.onchange();
    window.saveValues=async()=>{
      let blob; const create=URL.createObjectURL,click=HTMLAnchorElement.prototype.click;
      URL.createObjectURL=b=>{blob=b;return create.call(URL,b);};
      HTMLAnchorElement.prototype.click=function(){if(!this.download)click.call(this);};
      try{document.querySelector('#save-fit').click();}
      finally{URL.createObjectURL=create;HTMLAnchorElement.prototype.click=click;}
      return Object.assign({},...JSON.parse(await blob.text()).instrument.nodes.map(n=>n.parameters));
    };
    const before=await saveValues();
    const button=document.querySelector('.bloom-timing-button');
    const panel=document.querySelector('#bloom-timing-popup');
    const input=panel.querySelector('input');
    check(!panel.matches(':popover-open'),'Initially open');
    check(document.querySelectorAll('#bloom-controls fieldset').length===2,'Missing control groups');
    check(document.querySelector('#bloom-hold-controls [role=status]').textContent==='','Idle hold text clutter');
    button.scrollIntoView({block:'center'}); button.click();
    await new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)));
    check(panel.matches(':popover-open')&&button.getAttribute('aria-expanded')==='true','Did not open');
    check(document.activeElement===input,'Keyboard focus missing');
    input.value=.5;input.dispatchEvent(new Event('input',{bubbles:true}));
    const after=await saveValues();
    for(const key of Object.keys(before)) {
      const want=key==='bloom_rate'?before[key]*2**(-.5):key==='body_brightness'?before[key]-3:
        key==='body_excitation_centre'?before[key]*2**(-.125):before[key];
      check(Math.abs(after[key]-want)<1e-8,'Unexpected parameter change '+key);
    }
    check(document.querySelectorAll('.bloom-meta-target').length===3,'Missing real slider highlights');
    check(panel.querySelector('[data-preview-key=bloom_rate] output').textContent.includes('→'),'Missing live preview');
    panel.querySelector('[data-action=reset]').click();
    check(JSON.stringify(await saveValues())===JSON.stringify(before),'Reset did not restore baseline');
    input.value=.25;input.dispatchEvent(new Event('input',{bubbles:true}));
    const edited=await saveValues();panel.querySelector('[data-action=centre]').click();
    check(input.value==='0'&&JSON.stringify(await saveValues())===JSON.stringify(edited),'Recentre changed sound');
    input.value=.25;input.dispatchEvent(new Event('input',{bubbles:true}));
    const rect=panel.getBoundingClientRect();
    check(rect.left>=0&&rect.right<=innerWidth&&rect.top>=0&&rect.bottom<=innerHeight,'Popup outside viewport');
    return {parameterCount:Object.keys(before).length,preview:panel.querySelector('.bloom-timing-preview').textContent};
  })()`);
  await mkdir('build/bloom-timing-ui',{recursive:true});
  const shot = await call('Page.captureScreenshot',{format:'png'});
  await writeFile('build/bloom-timing-ui/popover.png',Buffer.from(shot.data,'base64'));
  await call('Input.dispatchKeyEvent',{type:'keyDown',key:'Escape',code:'Escape',windowsVirtualKeyCode:27});
  await call('Input.dispatchKeyEvent',{type:'keyUp',key:'Escape',code:'Escape',windowsVirtualKeyCode:27});
  assert.equal(await evaluate(`document.querySelector('#bloom-timing-popup').matches(':popover-open')`),false);
  assert.equal(await evaluate(`document.querySelectorAll('.bloom-meta-target').length`),0);
  await evaluate(`document.querySelector('.bloom-timing-button').click()`);
  await call('Input.dispatchMouseEvent',{type:'mousePressed',x:1450,y:80,button:'left',clickCount:1});
  await call('Input.dispatchMouseEvent',{type:'mouseReleased',x:1450,y:80,button:'left',clickCount:1});
  assert.equal(await evaluate(`document.querySelector('#bloom-timing-popup').matches(':popover-open')`),false);
  await evaluate(`document.querySelector('.bloom-timing-button').click()`);
  await call('Emulation.setDeviceMetricsOverride',{width:680,height:500,deviceScaleFactor:1,mobile:false});
  assert.equal(await evaluate(`new Promise(resolve=>requestAnimationFrame(()=>{
    const r=document.querySelector('#bloom-timing-popup').getBoundingClientRect();
    resolve(r.left>=0&&r.right<=innerWidth&&r.top>=0&&r.bottom<=innerHeight);
  }))`),true);
  await evaluate(`(async()=>{
    window.oldPopup=document.querySelector('#bloom-timing-popup');
    const selector=document.querySelector('#instrument-calibration');
    selector.value='crash-standard';await selector.onchange();
  })()`);
  assert.equal(await evaluate(`!oldPopup.isConnected&&document.querySelectorAll('#bloom-timing-popup').length===1&&
    document.querySelectorAll('.bloom-meta-target').length===0`),true);
  assert.equal(await evaluate(`(async()=>{
    const panel=document.querySelector('#bloom-timing-popup'),input=panel.querySelector('input');
    document.querySelector('.bloom-timing-button').click();
    input.dispatchEvent(new PointerEvent('pointerdown',{bubbles:true}));
    input.value=.1;input.dispatchEvent(new Event('input',{bubbles:true}));
    input.dispatchEvent(new Event('change',{bubbles:true}));
    await new Promise(r=>setTimeout(r,450));
    const status=document.querySelector('.decay-hold [role=status]');
    const ran=status.textContent.startsWith('Holding decay');
    document.querySelector('.decay-hold button').click();
    const cancelled=status.textContent.includes('Cancelled');
    const rate=document.querySelector('[data-fit-key=bloom_rate] input');
    rate.value=0;rate.dispatchEvent(new Event('input',{bubbles:true}));
    return ran&&cancelled&&input.disabled&&input.value==='0'&&
      panel.querySelector('.bloom-meta-status').textContent.includes('Enable Diffusion strength');
  })()`),true);
  assert.deepEqual(errors,[]);
  console.log('Bloom popover: preview, exact parameter changes, reset, recentre, dismissal, resize, rebuild, Hold decay and zero-strength state passed', result);
} finally {
  socket.close(); await fetch(endpoint+'/json/close/'+page.id);
}
