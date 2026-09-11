// Developer diagnostic: render a Plotly JSON figure in an isolated, silent tab.
import {readFile, writeFile} from "node:fs/promises";

const [input, output] = process.argv.slice(2);
if (!input || !output) throw Error("Usage: plotly_png.mjs input.plotly.json output.png");
const endpoint = "http://127.0.0.1:9223";
const page = await (await fetch(`${endpoint}/json/new?about:blank`, {method: "PUT"})).json();
const socket = new WebSocket(page.webSocketDebuggerUrl);
await new Promise(resolve => socket.addEventListener("open", resolve, {once: true}));
let id = 0;
const pending = new Map();
socket.onmessage = event => {
  const message = JSON.parse(event.data);
  if (message.id) {
    pending.get(message.id)?.(message);
    pending.delete(message.id);
  }
};
const call = (method, params = {}) => new Promise(resolve => {
  pending.set(++id, resolve);
  socket.send(JSON.stringify({id, method, params}));
});
try {
  const code = await readFile("workbench/web/node_modules/plotly.js-dist-min/plotly.min.js", "utf8");
  await call("Runtime.evaluate", {expression: code});
  const figure = JSON.parse(await readFile(input, "utf8"));
  const response = await call("Runtime.evaluate", {
    expression: `(async()=>{
      const f=${JSON.stringify(figure)};
      const div=document.body.appendChild(document.createElement('div'));
      await Plotly.newPlot(div,f.data,f.layout);
      return Plotly.toImage(div,{format:'png',width:f.layout.width,height:f.layout.height});
    })()`, awaitPromise: true, returnByValue: true,
  });
  if (response.error || response.result.exceptionDetails) throw Error(JSON.stringify(response));
  const encoded = response.result.result.value;
  await writeFile(output, Buffer.from(encoded.split(",")[1], "base64"));
} finally {
  socket.close();
  await fetch(`${endpoint}/json/close/${page.id}`);
}
