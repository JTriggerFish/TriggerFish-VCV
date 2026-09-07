// Render a local Plotly specification in a disposable debug tab, then close it.
import {readFile, writeFile} from "node:fs/promises";
import {resolve} from "node:path";
const directory = resolve(process.argv[2]);
const figure = JSON.parse(await readFile(resolve(directory,"inspection.plotly.json"),"utf8"));
const script = await readFile("build/workbench-wasm/site/vendor/plotly.min.js","utf8");
const endpoint = "http://127.0.0.1:9223";
const page = await (await fetch(endpoint+"/json/new?about:blank",{method:"PUT"})).json();
const socket = new WebSocket(page.webSocketDebuggerUrl);
await new Promise(resolve=>socket.addEventListener("open",resolve,{once:true}));
let sequence=0;
const pending=new Map();
socket.onmessage=event=>{const m=JSON.parse(event.data);if(pending.has(m.id)){
  const {resolve,reject}=pending.get(m.id);pending.delete(m.id);
  m.error?reject(Error(m.error.message)):resolve(m.result);
}};
const call=(method,params={})=>new Promise((resolve,reject)=>{
  const id=++sequence;pending.set(id,{resolve,reject});socket.send(JSON.stringify({id,method,params}));
});
try {
  await call("Emulation.setDeviceMetricsOverride",{width:1450,height:1150,deviceScaleFactor:1,mobile:false});
  const result=await call("Runtime.evaluate",{expression:
    `${script}\n document.body.style.margin='0'; document.body.innerHTML='<div id="plot"></div>'; Plotly.newPlot('plot',${JSON.stringify(figure.data)},${JSON.stringify(figure.layout)})`,awaitPromise:true});
  if(result.exceptionDetails)throw Error(JSON.stringify(result.exceptionDetails));
  await call("Runtime.evaluate", {expression:
    "document.fonts.ready.then(() => new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve))))",
    awaitPromise:true});
  const image=await call("Page.captureScreenshot",{format:"png"});
  await writeFile(resolve(directory,"inspection.png"),Buffer.from(image.data,"base64"));
} finally {
  socket.close();await fetch(endpoint+"/json/close/"+page.id);
}
