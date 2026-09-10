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


const click = async (x,y,count=1) => {
  await call('Input.dispatchMouseEvent',{type:'mouseMoved',x,y});
  for(let n=1;n<=count;n++){
    await call('Input.dispatchMouseEvent',{type:'mousePressed',x,y,button:'left',clickCount:n});
    await call('Input.dispatchMouseEvent',{type:'mouseReleased',x,y,button:'left',clickCount:n});
  }
};
const count = () => evaluate("document.querySelectorAll('#modal-editor .modal-node').length");
try {
  await call('Runtime.enable');await call('Page.enable');await call('Network.enable');
  await call('Network.setCacheDisabled',{cacheDisabled:true});
  await call('Emulation.setDeviceMetricsOverride',{width:1900,height:1200,deviceScaleFactor:1,mobile:false});
  await call('Page.addScriptToEvaluateOnNewDocument',{source:
    `Object.defineProperty(navigator,"requestMIDIAccess",{value:undefined});
    window.modalObservers=new Set();
    window.ResizeObserver=class extends ResizeObserver {
      observe(target,options){if(target.id==='modal-editor')modalObservers.add(this);super.observe(target,options);}
      disconnect(){modalObservers.delete(this);super.disconnect();}
    };`});
  await call('Page.navigate',{url:'http://127.0.0.1:8765/'});
  await evaluate(`new Promise((resolve,reject)=>{
    const deadline=performance.now()+45000;
    const poll=()=>document.querySelector('#instrument-calibration')?.options.length>2
      ?resolve(true):performance.now()>deadline?reject(Error('Not ready')):setTimeout(poll,100);poll();
  })`);
  await evaluate(`(async()=>{
    const s=document.querySelector('#instrument-calibration');s.value='gong-standard';await s.onchange();
    document.querySelector('#modal-clear').click();
    document.querySelector('#modal-editor').scrollIntoView({block:'center'});
    await new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)));
  })()`);
  const position=await evaluate(`(()=>{
    const r=document.querySelector('#modal-editor svg').getBoundingClientRect();
    return {x:r.left+r.width*.4,y:r.top+r.height*.4};
  })()`);
  await click(position.x,position.y,2);
  assert.equal(await count(),1,'Double-click empty plot should insert one mode');
  const handle=await evaluate(`(()=>{
    const r=document.querySelector('#modal-editor .modal-node').getBoundingClientRect();
    return {x:r.left+r.width/2,y:r.top+r.height/2};
  })()`);
  await click(handle.x,handle.y,2);
  assert.equal(await count(),0,'Real double-click on a repainted handle should remove it');
  const quiet=await evaluate(`(()=>{
    const r=document.querySelector('#modal-editor svg').getBoundingClientRect();
    return {x:r.left+r.width*.4,y:r.top+r.height*.8};
  })()`);
  await click(quiet.x,quiet.y,2);
  assert.ok(await evaluate(`parseFloat(document.querySelector('[data-fit-key=resolved_level_0] output').textContent)<-60`),
    'Quiet mode insertion must respect clicked level, not impose -48 dB');
  await click(quiet.x,quiet.y,2);
  assert.equal(await count(),0);
  await click(position.x,position.y,2);
  await call('Input.dispatchMouseEvent',{type:'mousePressed',...position,button:'left',clickCount:1});
  await call('Input.dispatchMouseEvent',{type:'mouseMoved',x:position.x+25,y:position.y-15,buttons:1});
  await call('Input.dispatchMouseEvent',{type:'mouseReleased',x:position.x+25,y:position.y-15,button:'left',clickCount:1});
  const moved=await evaluate(`(()=>{
    const r=document.querySelector('#modal-editor .modal-node').getBoundingClientRect();
    return {x:r.left+r.width/2,y:r.top+r.height/2};
  })()`);
  assert.ok(moved.x>position.x+15&&moved.y<position.y-5,'Centre dragging must still work');
  await call('Input.dispatchKeyEvent',{type:'keyDown',key:'Delete',code:'Delete',windowsVirtualKeyCode:46});
  await call('Input.dispatchKeyEvent',{type:'keyUp',key:'Delete',code:'Delete',windowsVirtualKeyCode:46});
  assert.equal(await count(),0,'Keyboard deletion should work after a pointer selection');
  await evaluate(`document.querySelector('#modal-tool-paint').click()`);
  await click(position.x,position.y,2);
  assert.equal(await count(),1,'Double-click in Paint mode should not delete the newly painted mode');
  await evaluate(`(async()=>{
    const s=document.querySelector('#instrument-calibration');s.value='gong-standard';await s.onchange();
    document.querySelector('#modal-editor').scrollIntoView({block:'center'});
    await new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)));
  })()`);
  assert.equal(await count(),32);
  assert.equal(await evaluate('modalObservers.size'),1,'Preset rebuild leaked a modal observer');
  const empty=await evaluate(`(()=>{
    const r=document.querySelector('#modal-editor svg').getBoundingClientRect();
    return {x:r.left+r.width*.15,y:r.top+r.height*.15};
  })()`);
  await click(empty.x,empty.y,2);
  assert.equal(await count(),32,'Full editor should not overwrite an existing mode');
  assert.equal(await evaluate(`!document.querySelector('#error-banner').hidden&&
    document.querySelector('#error-message').textContent.includes('All 32 modal handles')`),true,
    'Full capacity must not silently fail');
  await evaluate(`(async()=>{
    const s=document.querySelector('#instrument-calibration');s.value='kick-standard';await s.onchange();
  })()`);
  assert.equal(await evaluate('modalObservers.size'),0,'Recipe switch leaked a modal observer');
  assert.deepEqual(errors,[]);
  console.log('Real pointer interactions: insert, remove, drag, keyboard delete, painting and capacity notice pass');
} finally {
  socket.close();await fetch(endpoint+'/json/close/'+page.id);
}
