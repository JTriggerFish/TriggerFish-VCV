// Silent disposable browser profile; does not touch the user's playing tab.
import {writeFile} from 'node:fs/promises';
const endpoint = 'http://127.0.0.1:9223';
const page = await (await fetch(endpoint + '/json/new?about:blank', {method:'PUT'})).json();
const socket = new WebSocket(page.webSocketDebuggerUrl);
await new Promise(r => socket.addEventListener('open', r, {once:true}));
let id = 0;
const pending = new Map();
socket.onmessage = e => {
  const m = JSON.parse(e.data), p = pending.get(m.id);
  if (!p) return;
  pending.delete(m.id);
  m.error ? p.reject(Error(m.error.message)) : p.resolve(m.result);
};
const call = (method, params = {}) => new Promise((resolve, reject) => {
  pending.set(++id, {resolve, reject}); socket.send(JSON.stringify({id, method, params}));
});
const evaluate = async expression => {
  const result = await call('Runtime.evaluate', {expression, awaitPromise:true, returnByValue:true});
  if (result.exceptionDetails) throw Error(JSON.stringify(result.exceptionDetails));
  return result.result.value;
};
try {
  await call('Page.enable');
  await call('Page.addScriptToEvaluateOnNewDocument', {source:`
    Object.defineProperty(navigator,'requestMIDIAccess',{value:undefined});
    window.profileRows=[];
    window.Worker=class extends Worker {
      constructor(url, options) {
        super(url, options); this.url=String(url); this.starts=new Map();
        this.addEventListener('message',({data})=>{
          const key=data.kind+':'+data.generation, start=this.starts.get(key);
          if(start!==undefined){profileRows.push({worker:this.url,kind:data.kind,
            wallMs:performance.now()-start,dspMs:data.elapsedMs,
            received:performance.now(),preview:data.preview});if(!data.preview)this.starts.delete(key);}
        });
      }
      postMessage(data,...rest){this.starts.set(data.kind+':'+data.generation,performance.now());super.postMessage(data,...rest);}
    };`});
  await call('Page.navigate', {url:'http://127.0.0.1:8765/'});
  const results = await evaluate(`(async()=>{
    const wait=async test=>{const until=performance.now()+60000;while(!test()){
      if(performance.now()>until)throw Error('Profile timed out');await new Promise(r=>setTimeout(r,20));}};
    await wait(()=>document.querySelector('#instrument-calibration')?.options.length>2);
    const select=document.querySelector('#instrument-calibration'); select.value='gong-standard';await select.onchange();
    await wait(()=>document.querySelector('#status').textContent==='Ready');
    const {SpectrogramView}=await import('./spectrogram.mjs');
    const draw=SpectrogramView.prototype.draw;
    SpectrogramView.prototype.draw=function(...args){const start=performance.now();
      const result=draw.apply(this,args);profileRows.push({drawMs:performance.now()-start});
      if(this.synthesis?.incomplete){
        if(this.synthesis.frames<=this.synthesis.writeFrames)throw Error('Previous spectrogram tail was lost');
        if(!window.previewPng && this.synthesis.writeFrames>100)
          window.previewPng=this.canvas.toDataURL('image/png');
      }
      return result;};
    const input=document.querySelector('[data-fit-key=field_motion_depth] input[type=range]');
    const results=[];
    for(const edits of [1,20]){
      profileRows.length=0;const start=performance.now();
      for(let i=0;i<edits;i++){input.value=1.5+.1*(i%3);input.dispatchEvent(new Event('input',{bubbles:true}));
        if(edits>1)await new Promise(r=>setTimeout(r,40));}
      const endEdit=performance.now();await wait(()=>document.querySelector('#status').textContent==='Ready');
      const first=profileRows.find(r=>r.worker==='analysis_worker.mjs'&&r.kind==='synthesis');
      results.push({edits,totalMs:performance.now()-start,afterLastEditMs:performance.now()-endEdit,
        firstAnalysisMs:first?.received-start,rows:profileRows.filter(r=>r.worker!=='preparation_worker.mjs')});
    }
    return results;
  })()`);
  console.log(JSON.stringify(results.map(({rows, ...timing}) => ({...timing,
    fullRenderMs: rows.find(r => r.worker === 'render_worker.mjs' && !r.preview)?.dspMs,
    maximumAnalysisMs: Math.max(...rows.filter(r => r.worker === 'analysis_worker.mjs').map(r => r.wallMs)),
    updates: rows.filter(r => r.drawMs !== undefined).length,
  })), null, 2));
  if (process.argv[2]) {
    const png = await evaluate(`(()=>{
      return (window.previewPng ?? document.querySelector('#spectrogram').toDataURL('image/png')).split(',')[1];
    })()`);
    await writeFile(process.argv[2], Buffer.from(png, 'base64'));
  }
} finally {
  socket.close(); await fetch(endpoint+'/json/close/'+page.id);
}
