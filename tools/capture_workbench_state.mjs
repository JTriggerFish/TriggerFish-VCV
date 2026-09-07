// Read the open workbench's module state without clicking, reloading or editing.
import {writeFile} from "node:fs/promises";
const tabs=await(await fetch("http://127.0.0.1:9223/json/list")).json();
const tab=tabs.find(t=>t.type==="page" && /:8765\/$/.test(t.url));
if(!tab)throw Error("No open workbench tab");
const socket=new WebSocket(tab.webSocketDebuggerUrl);
await new Promise(resolve=>socket.addEventListener("open",resolve,{once:true}));
let id=0;const pending=new Map();
socket.onmessage=event=>{const r=JSON.parse(event.data),p=pending.get(r.id);if(p){pending.delete(r.id);r.error?p.reject(Error(r.error.message)):p.resolve(r.result);}};
const call=(method,params={})=>new Promise((resolve,reject)=>{pending.set(++id,{resolve,reject});socket.send(JSON.stringify({id,method,params}));});
try{
  const fn=await call("Runtime.evaluate",{expression:"document.getElementById('save-fit').onclick"});
  const props=await call("Runtime.getProperties",{objectId:fn.result.objectId});
  const scopes=props.internalProperties.find(p=>p.name==="[[Scopes]]").value;
  const list=await call("Runtime.getProperties",{objectId:scopes.objectId});
  let state,engine;
  for(const scope of list.result.filter(p=>/^\d+$/.test(p.name))){
    const values=await call("Runtime.getProperties",{objectId:scope.value.objectId});
    state??=values.result.find(p=>p.name==="state")?.value;
    engine??=values.result.find(p=>p.name==="engine")?.value;
  }
  if(!state?.objectId||!engine?.objectId)throw Error("Cannot inspect workbench state");
  const captured=await call("Runtime.callFunctionOn",{objectId:state.objectId,
    functionDeclaration:`async function(engine){const {snapshotState}=await import('${tab.url}state.mjs');return snapshotState(this,'Before fixed-strike kick cleanup',engine.macros);}`,
    arguments:[{objectId:engine.objectId}],awaitPromise:true,returnByValue:true});
  if(captured.exceptionDetails)throw Error(JSON.stringify(captured.exceptionDetails));
  await writeFile(process.argv[2],JSON.stringify(captured.result.value,null,2)+"\n");
  const shot=await call("Page.captureScreenshot",{format:"png",captureBeyondViewport:true});
  await writeFile(process.argv[2]+".png",Buffer.from(shot.data,"base64"));
  console.log(captured.result.value.renderer,captured.result.value.controls.event);
}finally{socket.close();}
