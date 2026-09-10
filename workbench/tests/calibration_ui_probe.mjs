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
  await call("Page.navigate",{url:"http://127.0.0.1:8765/"});
  await evaluate(`new Promise((resolve,reject)=>{
    const deadline=performance.now()+45000;
    const poll=()=>{
      if(document.getElementById('instrument-calibration')?.options.length>2)resolve(true);
      else if(performance.now()>deadline)reject(Error('Targets did not load'));
      else setTimeout(poll,100);
    };poll();
  })`);
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
                                 ['bloom_energy_sensitivity','Energy sensitivity']]) {
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
  if(errors.length)throw Error(errors.join('\n'));
} finally {
  socket.close();await fetch(endpoint+"/json/close/"+page.id);
}
