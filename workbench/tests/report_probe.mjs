// Inspect the local fitting report in a temporary tab, not the user's workbench.
import { writeFile } from "node:fs/promises";

const endpoint = "http://127.0.0.1:9223";
const report = process.env.TF_FIT_REPORT ?? "gong-review";
if (!/^[a-z-]+(?:\/[a-z-]+)*$/.test(report)) throw new Error("Invalid report name");
const url = `http://127.0.0.1:8765/${report}/`;
const page = await (await fetch(`${endpoint}/json/new?about:blank`, { method: "PUT" })).json();
const socket = new WebSocket(page.webSocketDebuggerUrl);
await new Promise((resolve, reject) => {
  socket.addEventListener("open", resolve, { once: true });
  socket.addEventListener("error", reject, { once: true });
});
let nextId = 1;
let loaded = false;
const pending = new Map();
socket.onmessage = event => {
  const response = JSON.parse(event.data);
  if (response.method === "Page.loadEventFired") loaded = true;
  if (!pending.has(response.id)) return;
  const [resolve, reject] = pending.get(response.id);
  pending.delete(response.id);
  if (response.error) reject(new Error(response.error.message));
  else resolve(response.result);
};
const call = (method, params = {}) => new Promise((resolve, reject) => {
  const id = nextId++;
  pending.set(id, [resolve, reject]);
  socket.send(JSON.stringify({ id, method, params }));
});
try {
  await call("Emulation.setDeviceMetricsOverride", {
    width: 1500, height: report.startsWith("kick-") ? 2500 : 1750, deviceScaleFactor: 1, mobile: false,
  });
  await call("Runtime.enable");
  await call("Page.enable");
  await call("Page.navigate", { url });
  const deadline = Date.now() + 20000;
  while (!loaded) {
    if (Date.now() > deadline) throw new Error("Report navigation timed out");
    await new Promise(resolve => setTimeout(resolve, 50));
  }
  const ready = await call("Runtime.evaluate", {
    expression: `new Promise((resolve,reject)=>{
      const start=performance.now();
      const poll=()=>{
        const plots=[...document.querySelectorAll('.plotly-graph-div')];
        if(plots.length>=2&&plots.every(p=>p._fullLayout))resolve(true);
        else if(performance.now()-start>20000)reject(Error('report did not load'));
        else setTimeout(poll,100);
      };poll();
    })`, awaitPromise: true, returnByValue: true,
  });
  if (ready.result.value !== true) throw new Error(JSON.stringify(ready));
  const checked = await call("Runtime.evaluate", {
    expression: `(async()=>{
      const fit=await(await fetch(document.querySelector('a[download]').getAttribute('href'))).json();
      const audio=await Promise.all(['reference','candidate'].map(async name=>
        (await(await fetch(document.querySelector('[data-audio="'+name+'"]').dataset.file??name+'.wav')).arrayBuffer()).byteLength>100000));
      const active=fit.instrument.nodes.flatMap(n=>Object.entries(n.parameters))
        .filter(([k,v])=>/^body_decay_active_[1-6]$/.test(k)&&v>=.5);
      document.querySelector('[data-audio="candidate"]').click();
      await new Promise(resolve=>setTimeout(resolve,1500));
      const protectedPlayback=document.getElementById('audio-status').textContent.includes('limiter enabled');
      document.getElementById('stop').click();
      return {audio:audio.every(Boolean),twoPointDecay:active.length===0,protectedPlayback,
        status:document.getElementById('audio-status').textContent};
    })()`, awaitPromise: true, returnByValue: true, userGesture: true,
  });
  console.log(JSON.stringify(checked.result.value));
  if (!checked.result.value || !["audio", "twoPointDecay", "protectedPlayback"].every(key => checked.result.value[key]))
    throw new Error("Report checks failed");
  const screenshot = await call("Page.captureScreenshot", { format: "png", captureBeyondViewport: true });
  await writeFile(`build/${report.replaceAll("/", "-")}.png`, Buffer.from(screenshot.data, "base64"));
} finally {
  await new Promise(resolve => {
    socket.addEventListener("close", resolve, { once: true });
    socket.close();
  });
  await fetch(`${endpoint}/json/close/${page.id}`);
}
