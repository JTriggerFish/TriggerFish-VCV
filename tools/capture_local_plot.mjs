// Capture a self-contained local diagnostic in a temporary tab, then close it.
import {resolve} from "node:path";
import {pathToFileURL} from "node:url";
import {writeFile} from "node:fs/promises";
const endpoint="http://127.0.0.1:9223";
const page=await(await fetch(`${endpoint}/json/new?about:blank`,{method:"PUT"})).json();
const socket=new WebSocket(page.webSocketDebuggerUrl);
await new Promise(resolve=>socket.addEventListener("open",resolve,{once:true}));
let id=0;
const pending=new Map();
socket.onmessage=event=>{
  const response=JSON.parse(event.data), waiter=pending.get(response.id);
  if(waiter){pending.delete(response.id); response.error?waiter.reject(Error(response.error.message)):waiter.resolve(response.result);}
};
const call=(method,params={})=>new Promise((resolve,reject)=>{
  pending.set(++id,{resolve,reject});socket.send(JSON.stringify({id,method,params}));
});
try{
  await call("Emulation.setDeviceMetricsOverride",{width:1520,height:1020,deviceScaleFactor:1,mobile:false});
  await call("Page.navigate",{url:pathToFileURL(resolve(process.argv[2])).href});
  const deadline=Date.now()+20000;
  for(;;){
    const ready=await call("Runtime.evaluate",{expression:"!!document.querySelector('.plotly-graph-div')?._fullLayout",returnByValue:true});
    if(ready.result?.value)break;
    if(Date.now()>deadline)throw Error("Plot did not render");
    await new Promise(resolve=>setTimeout(resolve,100));
  }
  const shot=await call("Page.captureScreenshot",{format:"png",captureBeyondViewport:true});
  await writeFile(resolve(process.argv[3]),Buffer.from(shot.data,"base64"));
}finally{
  socket.close();await fetch(`${endpoint}/json/close/${page.id}`);
}
