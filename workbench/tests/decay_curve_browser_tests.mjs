// Optional developer CDP check in a disposable silent tab, not the user's page.
import { readFile } from "node:fs/promises";
const endpoint = process.env.TF_CDP_URL ?? "http://127.0.0.1:9223";
const page = await (await fetch(`${endpoint}/json/new?about:blank`, { method: "PUT" })).json();
const socket = new WebSocket(page.webSocketDebuggerUrl);
await new Promise(resolve => socket.addEventListener("open", resolve, { once: true }));
let sequence = 0;
const pending = new Map();
socket.onmessage = event => {
  const message = JSON.parse(event.data);
  const item = pending.get(message.id);
  if (!item) return;
  pending.delete(message.id);
  message.error ? item.reject(Error(message.error.message)) : item.resolve(message.result);
};
const call = (method, params = {}) => new Promise((resolve, reject) => {
  const id = ++sequence; pending.set(id, { resolve, reject });
  socket.send(JSON.stringify({ id, method, params }));
});
try {
  const geometry = (await readFile(new URL("../web/decay_curve_geometry.mjs", import.meta.url), "utf8")).replaceAll("export ", "");
  const editor = (await readFile(new URL("../web/decay_curve_editor.mjs", import.meta.url), "utf8"))
    .replace(/^import .*\n/, "").replace("export class", "class");
  const result = await call("Runtime.evaluate", { awaitPromise: true, returnByValue: true,
    expression: `(async () => {
      ${geometry}\n${editor}
      document.body.innerHTML = '<style>svg{width:100%;height:220px}</style><div id="editor" style="width:240px"></div>';
      let points = [{slot:0,x:40,y:Math.log2(29),fixed:true},{slot:7,x:15000,y:0,fixed:true}];
      const initial = points.map(p => ({...p}));
      const control = new DecayCurveEditor(document.getElementById('editor'), {
        minimumFrequency:40,maximumFrequency:15000, minimumLogSeconds:Math.log2(.02),maximumLogSeconds:Math.log2(30),
        yTicks:[],points:()=>points, select:()=>{},
        setPoint:(slot,x,y)=>{ points=points.map(p=>p.slot===slot?{...p,x,y}:p); },
        replace:p=>{points=p;}, insert:()=>null, remove:()=>{},reset:()=>{}
      });
      await new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)));
      const matrix = control.svg.getScreenCTM();
      const origin = {x:control.xPosition(40),y:control.yPosition(initial[0].y)};
      const client = {clientX:matrix.a*origin.x+matrix.c*origin.y+matrix.e,
        clientY:matrix.b*origin.x+matrix.d*origin.y+matrix.f};
      const mapped = control.eventPosition(client);
      if (Math.abs(mapped.y-origin.y)>1e-8) throw Error('SVG pointer coordinates differ');
      if (origin.y<18 || origin.y>186) throw Error('29-second endpoint outside graph');
      control.drag={slot:0,start:origin};
      control.dragPoint(origin);
      if (points[0].y!==initial[0].y) throw Error('No-motion drag changes value');
      control.dragPoint({...origin,y:origin.y+10});
      const normal = initial[0].y-points[0].y;
      points=initial.map(p=>({...p}));
      control.dragPoint({...origin,y:origin.y+10},true);
      const fine = initial[0].y-points[0].y;
      if (Math.abs(normal/fine-10)>1e-8) throw Error('Fine drag ratio differs');
      document.getElementById('editor').style.width='900px';
      await new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)));
      if (Math.abs(control.width-900)>1) throw Error('Curve does not use resized width');
      control.destroy();
      return {coordinateMapping:true,noJump:true,fineRatio:normal/fine,fullWidth:control.width};
    })()` });
  if (result.exceptionDetails) throw Error(JSON.stringify(result.exceptionDetails));
  console.log(JSON.stringify(result.result.value));
} finally {
  socket.close();
  await fetch(`${endpoint}/json/close/${page.id}`);
}
