// Integration check in a temporary tab; never reload the user's workbench.
import { writeFile } from "node:fs/promises";
const endpoint = "http://127.0.0.1:9223";
const page = await (await fetch(endpoint + "/json/new?about:blank", {method:"PUT"})).json();
const socket = new WebSocket(page.webSocketDebuggerUrl);
await new Promise(resolve => socket.addEventListener("open", resolve, {once:true}));
let nextId = 0;
const pending = new Map(), errors = [];
socket.onmessage = event => {
  const message = JSON.parse(event.data);
  if (message.method === "Runtime.exceptionThrown")
    errors.push(message.params.exceptionDetails.text);
  const waiter = pending.get(message.id);
  if (waiter) { pending.delete(message.id); waiter(message.result); }
};
const call = (method, params = {}) => new Promise(resolve => {
  const id = ++nextId; pending.set(id, resolve);
  socket.send(JSON.stringify({id, method, params}));
});
const evaluate = async expression => {
  const result = await call("Runtime.evaluate", {
    expression, awaitPromise:true, returnByValue:true, userGesture:true,
  });
  if (result.exceptionDetails) throw Error(JSON.stringify(result.exceptionDetails));
  return result.result.value;
};
try {
  await call("Runtime.enable");
  await call("Page.enable");
  await call("Page.navigate", {url:"http://127.0.0.1:8765/"});
  await evaluate(`new Promise((resolve,reject)=>{
    const deadline=performance.now()+30000;
    const poll=()=>{
      const target=document.getElementById('instrument-calibration');
      if([...target?.options??[]].some(o=>o.value==='kick-standard')) resolve(true);
      else if(performance.now()>deadline) reject(Error('Targets did not load'));
      else setTimeout(poll,100);
    };poll();
  })`);
  await evaluate(`(async()=>{
    const target=document.getElementById('instrument-calibration');
    target.value='kick-standard'; await target.onchange();
    return true;
  })()`);
  const result = await evaluate(`(async()=>{
    const {PercussionEngine}=await import('./engine.mjs');
    const engine=await PercussionEngine.create(44100);
    engine.setRecipe(engine.recipes.find(r=>r.key==='drum.kick.v1').index);
    const keys=engine.parameters.map(d=>d.key).sort();
    const rows=[...document.querySelectorAll('[data-kick-key]')];
    const actual=rows.map(r=>r.dataset.kickKey).sort();
    const base=keys.filter(k=>!/^resonance_(frequency|level|centre|edge)_[0-9]+$/.test(k));
    if(!base.every(k=>actual.includes(k)) ||
       keys.filter(k=>/^resonance_(frequency|level|centre|edge)_[0-9]+$/.test(k)).length!==64)
      throw Error('UI/DSP control mismatch');
    const generator=document.querySelector('#kick-modal-templates button');
    generator.click();
    const bars=document.querySelectorAll('#kick-modal-editor .modal-bar');
    if(bars.length!==16) throw Error('Membrane meta tool did not create editable modes');
    document.getElementById('kick-mode-clear').click();
    if(document.querySelectorAll('#kick-modal-editor .modal-bar').length) throw Error('Clear did not remove modes');
    generator.click();
    const eq=document.querySelector('[data-kick-key="equalizer_mode"] select');
    if(eq.value!=='0') throw Error('Kick reference start must bypass EQ');
    for(const mode of [2,1,0]){
      eq.value=mode; eq.dispatchEvent(new Event('change'));
      if(document.getElementById('kick-radiation-controls').hidden!==(mode===0) ||
         document.querySelector('[data-kick-key="colour_gain_db"]').hidden!==(mode!==1) ||
         document.getElementById('kick-multiband-controls').hidden!==(mode!==2))
        throw Error('EQ section visibility mismatch');
    }
    const family=document.getElementById('kick-implement');
    if(family.value!=='0.5') throw Error('Reference implement was not restored');
    document.getElementById('play-synthesis').click();
    await new Promise(resolve=>setTimeout(resolve,2500));
    engine.destroy();
    return {controls:keys.length, family:family.value,
      recipe:document.getElementById('instrument-recipe').selectedOptions[0].text,
      status:document.getElementById('status').textContent,
      live:document.getElementById('live-commit').textContent};
  })()`);
  if (errors.length) throw Error(errors.join("\n"));
  if (/error|failed|unavailable/i.test(result.status)) throw Error(result.status);
  await evaluate(`(async()=>{
    const {PercussionEngine}=await import('./engine.mjs');
    const {FitControls}=await import('./fit_controls.mjs');
    const engine=await PercussionEngine.create(44100,0);
    const controls=new FitControls({descriptors:engine.parameters,
      state:{macros:engine.parameters.map(d=>d.defaultValue)},onChange:()=>{}});
    document.getElementById('modal-editor').replaceChildren();
    controls.buildResolvedEditor();
    controls.resolvedEditor.select(0);
    document.querySelector('#modal-templates button').click();
    if(controls.resolvedEditor.selected!==null ||
       !document.querySelector('#modal-selection b').textContent.startsWith('No '))
      throw Error('Modal generation left a stale selection');
    const before=controls.state.macros.slice();
    controls.resolvedEditor.svg.dispatchEvent(new KeyboardEvent('keydown',{key:'Delete'}));
    if(JSON.stringify(before)!==JSON.stringify(controls.state.macros))
      throw Error('Delete removed an unselected generated mode');
    controls.resolvedEditor.resizeObserver.disconnect(); engine.destroy();
  })()`);
  console.log(JSON.stringify(result));
  const shot = await call("Page.captureScreenshot", {format:"png", captureBeyondViewport:true});
  await writeFile("build/kick-workbench.png", Buffer.from(shot.data, "base64"));
} finally {
  socket.close();
  await fetch(endpoint + "/json/close/" + page.id);
}
