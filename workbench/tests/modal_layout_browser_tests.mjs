// Silent disposable browser tab: validate layout, not audio or user edits.
import {writeFile} from "node:fs/promises";
import {checkDecayInteractions} from "./decay_editor_browser_checks.mjs";
const endpoint = "http://127.0.0.1:9223";
const page = await (await fetch(endpoint + "/json/new?about:blank", {method:"PUT"})).json();
const socket = new WebSocket(page.webSocketDebuggerUrl);
await new Promise(resolve => socket.addEventListener("open", resolve, {once:true}));
let next = 0;
const pending = new Map();
socket.onmessage = ({data}) => {
  const message = JSON.parse(data), task = pending.get(message.id);
  if (!task) return;
  pending.delete(message.id);
  message.error ? task.reject(Error(message.error.message)) : task.resolve(message.result);
};
const call = (method, params = {}) => new Promise((resolve, reject) => {
  pending.set(++next, {resolve, reject}); socket.send(JSON.stringify({id:next, method, params}));
});
const evaluate = async expression => {
  const result = await call("Runtime.evaluate", {expression, awaitPromise:true, returnByValue:true});
  if (result.exceptionDetails) throw Error(JSON.stringify(result.exceptionDetails));
  return result.result.value;
};
try {
  await call("Page.enable");
  await call("Page.addScriptToEvaluateOnNewDocument", {source:
    'Object.defineProperty(navigator,"requestMIDIAccess",{value:undefined});'});
  await call("Emulation.setDeviceMetricsOverride", {width:2560,height:1440,deviceScaleFactor:1,mobile:false});
  await call("Page.navigate", {url:"http://127.0.0.1:8765/"});
  await evaluate(`(async()=>{
    const until=performance.now()+30000;
    while(typeof document.querySelector('#instrument-calibration')?.onchange!=='function'){
      if(performance.now()>until)throw Error('Initialization timed out');
      await new Promise(r=>setTimeout(r,50));
    }
    const target=document.getElementById('instrument-calibration');
    target.value='gong-standard';await target.onchange();
    const ids=[...document.querySelectorAll('[id]')].map(e=>e.id);
    if(new Set(ids).size!==ids.length)throw Error('Duplicate element IDs');
    for(const id of ['modal-editor','modal-templates','strike-pad','spectrogram'])
      if(!document.getElementById(id).closest('.analysis'))throw Error('Not in analysis panel: '+id);
    for(const id of ['bloom-controls','field-drive-controls','field-tuning-controls','field-turbulence-controls',
      'field-beating-controls','field-drift-controls','field-motion-controls','field-blur-controls'])
      if(!document.getElementById(id).closest('#resonance-column'))throw Error('Not in right control column: '+id);
    for(const id of ['impact-controls','decay-editor','output-eq-editor'])
      if(!document.getElementById(id).closest('#excitation-column'))throw Error('Not in left control column: '+id);
    const routing=document.getElementById('routing-overview');
    if(routing.tagName!=='DETAILS' || routing.open)throw Error('Routing must start collapsed');
    routing.querySelector('summary').click();
    if(!routing.open)throw Error('Routing accordion did not open');
    routing.querySelector('summary').click();
    const input=key=>document.querySelector('[data-fit-key='+key+'] input');
    if(!input('field_wander_rate').disabled || input('field_motion_rate').disabled ||
       !input('field_phase_tilt').disabled)throw Error('Incorrect initial dependent states');
    input('field_wander_hz').value=.2;input('field_wander_hz').dispatchEvent(new Event('input'));
    if(input('field_wander_rate').disabled)throw Error('Drift speed did not enable');
    input('field_motion_depth').value=0;input('field_motion_depth').dispatchEvent(new Event('input'));
    if(!input('field_motion_rate').disabled || !input('field_motion_sharing').disabled)
      throw Error('Shimmer dependents did not disable');
    await target.onchange();
  })()`);
  for (const width of [2560, 1920, 1440]) {
    await call("Emulation.setDeviceMetricsOverride", {width,height:1100,deviceScaleFactor:1,mobile:false});
    console.log(await evaluate(`(async()=>{
      await new Promise(requestAnimationFrame);await new Promise(requestAnimationFrame);
      const left=document.querySelector('aside'),right=document.querySelector('.analysis');
      const columns=[...left.querySelectorAll('.control-column')];
      if(columns.length!==2)throw Error('Expected two control columns in aside');
      const [a,b]=columns.map(e=>e.getBoundingClientRect());
      if(b.left<a.right || Math.abs(a.top-b.top)>1)throw Error('Control columns are not side by side');
      if(right.scrollWidth>right.clientWidth+2)throw Error('Right column overflows at '+innerWidth+' '+JSON.stringify(
        [...right.querySelectorAll('*')].filter(e=>e.getClientRects().length && e.getBoundingClientRect().right>right.getBoundingClientRect().right)
          .slice(0,12).map(e=>[e.tagName,e.id,e.className?.baseVal??e.className,e.getBoundingClientRect().width])));
      if(left.scrollWidth>left.clientWidth+2)throw Error('Left column overflows at '+innerWidth);
      const before=left.scrollTop;right.scrollTop=600;
      if(left.scrollTop!==before || right.scrollTop!==600)throw Error('Column scrolling is not independent');
      return {width:innerWidth,left:left.clientWidth,right:right.clientWidth};
    })()`));
  }
  await call("Emulation.setDeviceMetricsOverride", {width:1920,height:1200,deviceScaleFactor:1,mobile:false});
  await evaluate(`(async()=>{const until=performance.now()+30000;
    while(document.getElementById('status').textContent!=='Ready'){
      if(performance.now()>until)throw Error('Render did not complete');
      await new Promise(r=>setTimeout(r,50));
    }})()`);
  await checkDecayInteractions(call, evaluate);
  for (const [name, expression] of [
    ["overview", "document.querySelector('.analysis').scrollTop=0"],
    ["controls", "document.getElementById('resonance-texture').scrollIntoView({block:'start'})"],
  ]) {
    await evaluate(`(async()=>{${expression};await new Promise(requestAnimationFrame);})()`);
    const png = await call("Page.captureScreenshot", {format:"png"});
    await writeFile(`build/modal-layout-${name}.png`, Buffer.from(png.data,"base64"));
  }
} finally {
  socket.close();await fetch(endpoint+"/json/close/"+page.id);
}
