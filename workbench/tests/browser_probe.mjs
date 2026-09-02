import { writeFile } from "node:fs/promises";

const arguments_ = process.argv.slice(2);
const endpoint = arguments_.find(argument => /^https?:\/\//.test(argument)) ??
  "http://127.0.0.1:9223";
const screenshot = arguments_.find(argument =>
  argument !== endpoint && !argument.startsWith("--"));
const testAudio = arguments_.includes("--trigger");
const profileUi = arguments_.includes("--profile");
const testControls = arguments_.includes("--controls");

async function pageTarget() {
  const targets = await (await fetch(`${endpoint}/json`)).json();
  return targets.find(target => target.type === "page");
}

function connect(url) {
  return new Promise((resolve, reject) => {
    const socket = new WebSocket(url);
    socket.addEventListener("open", () => resolve(socket), { once: true });
    socket.addEventListener("error", reject, { once: true });
  });
}

function protocol(socket) {
  let nextId = 1;
  const pending = new Map();
  socket.addEventListener("message", event => {
    const message = JSON.parse(event.data);
    if (!message.id || !pending.has(message.id)) return;
    const { resolve, reject } = pending.get(message.id);
    pending.delete(message.id);
    if (message.error) reject(new Error(message.error.message));
    else resolve(message.result);
  });
  return (method, params = {}) => new Promise((resolve, reject) => {
    const id = nextId++;
    pending.set(id, { resolve, reject });
    socket.send(JSON.stringify({ id, method, params }));
  });
}

const target = await pageTarget();
if (!target) throw new Error("no browser page target");
const socket = await connect(target.webSocketDebuggerUrl);
const call = protocol(socket);
await call("Runtime.enable");
const expression = `new Promise(resolve => {
  const deadline = performance.now() + 15000;
  const poll = () => {
    const status = document.getElementById("status")?.textContent ?? "";
    const canvas = document.getElementById("spectrogram");
    const pixels = canvas?.getContext("2d").getImageData(
      0, 0, canvas.width, canvas.height).data;
    let lit = 0;
    if (pixels) for (let i = 0; i < pixels.length; i += 64) {
      if (pixels[i] + pixels[i + 1] + pixels[i + 2] > 12) ++lit;
    }
    if (status === "Ready" && lit > 0 || performance.now() >= deadline) {
      resolve({ status, lit, width: canvas?.width, height: canvas?.height });
    } else setTimeout(poll, 100);
  };
  poll();
})`;
const evaluated = await call("Runtime.evaluate", {
  expression, awaitPromise: true, returnByValue: true,
});
const result = evaluated.result.value;
if (profileUi) {
  const profile = await call("Runtime.evaluate", {
    expression: `(() => {
      const range = document.getElementById("colour-range");
      const mode = document.getElementById("view-mode");
      const measure = action => {
        const started = performance.now(); action();
        return performance.now() - started;
      };
      const rangeDrawMs = [];
      for (let repeat = 0; repeat < 12; ++repeat) {
        rangeDrawMs.push(measure(() => {
          range.value = repeat % 2 ? 88 : 92;
          range.dispatchEvent(new Event("input"));
        }));
      }
      const modeDrawMs = [];
      for (const value of ["reference", "synthesis", "mirror", "difference",
                           "mirror", "vertical", "mirror"]) {
        modeDrawMs.push(measure(() => {
          mode.value = value; mode.dispatchEvent(new Event("change"));
        }));
      }
      const median = values => [...values].sort((a, b) => a - b)[
        Math.floor(values.length / 2)];
      return {
        pixels: document.getElementById("spectrogram").width *
          document.getElementById("spectrogram").height,
        rangeDrawMedianMs: median(rangeDrawMs),
        rangeDrawMaximumMs: Math.max(...rangeDrawMs),
        modeDrawMedianMs: median(modeDrawMs),
        modeDrawMaximumMs: Math.max(...modeDrawMs),
      };
    })()`,
    returnByValue: true,
  });
  result.profile = profile.result.value;
  const interaction = await call("Runtime.evaluate", {
    expression: `new Promise(resolve => {
      const longTasks = [];
      const observer = typeof PerformanceObserver === "undefined" ? null :
        new PerformanceObserver(list => {
          for (const entry of list.getEntries()) longTasks.push(entry.duration);
        });
      observer?.observe({ type: "longtask" });
      let maximumTimerGapMs = 0;
      let previous = performance.now();
      const timer = setInterval(() => {
        const now = performance.now();
        maximumTimerGapMs = Math.max(maximumTimerGapMs, now - previous);
        previous = now;
      }, 10);
      const slider = document.querySelector("#model-level input");
      const started = performance.now();
      slider.value = Number(slider.value) + .1;
      slider.dispatchEvent(new Event("input"));
      const poll = () => {
        if (document.getElementById("status").textContent === "Ready" &&
            performance.now() - started > 300) {
          clearInterval(timer); observer?.disconnect();
          resolve({
            elapsedMs: performance.now() - started,
            dspLabel: document.getElementById("render-time").textContent,
            maximumTimerGapMs,
            longTaskCount: longTasks.length,
            longestTaskMs: Math.max(0, ...longTasks),
          });
        } else setTimeout(poll, 20);
      };
      poll();
    })`,
    awaitPromise: true, returnByValue: true,
  });
  result.profile.interaction = interaction.result.value;
}
if (testAudio) {
  const button = await call("Runtime.evaluate", {
    expression: `(() => {
      const bounds = document.getElementById("play-synthesis").getBoundingClientRect();
      return { x: bounds.x + bounds.width / 2, y: bounds.y + bounds.height / 2 };
    })()`,
    returnByValue: true,
  });
  for (let strike = 0; strike < 3; ++strike) {
    await call("Input.dispatchMouseEvent", {
      type: "mousePressed", x: button.result.value.x, y: button.result.value.y,
      button: "left", clickCount: 1,
    });
    await call("Input.dispatchMouseEvent", {
      type: "mouseReleased", x: button.result.value.x, y: button.result.value.y,
      button: "left", clickCount: 1,
    });
    await new Promise(resolve => setTimeout(resolve, 120));
  }
  const audioProbe = await call("Runtime.evaluate", {
    expression: `new Promise(resolve => {
      const deadline = performance.now() + 15000;
      const poll = () => {
        const status = document.getElementById("status")?.textContent ?? "";
        const meter = document.getElementById("limiter")?.textContent ?? "";
        if (meter.includes("running") && /hits [1-9]\d*/.test(meter) ||
            status.includes("Error") || performance.now() >= deadline) {
          resolve({ status, meter });
        } else setTimeout(poll, 100);
      };
      poll();
    })`,
    awaitPromise: true, returnByValue: true,
  });
  const meterProbe = await call("Runtime.evaluate", {
    expression: `new Promise(resolve => setTimeout(() => resolve({
      status: document.getElementById("status")?.textContent ?? "",
      meter: document.getElementById("limiter")?.textContent ?? "",
    }), 700))`,
    awaitPromise: true, returnByValue: true,
  });
  result.audio = { ...audioProbe.result.value, afterStrikes: meterProbe.result.value };
}
if (testControls) {
  const controls = await call("Runtime.evaluate", {
    expression: `new Promise(resolve => {
      const master = document.getElementById("master");
      const hardness = document.getElementById("hardness");
      const colour = document.getElementById("colour-range");
      const mode = document.getElementById("view-mode");
      const model = document.querySelector("#model-level input");
      const decay = document.querySelector("#decay-controls input");
      const initial = {
        hardness: hardness.value, model: Number(model.value), decay: decay.value,
      };
      const reset = (element, changed, eventName = "input") => {
        element.value = changed;
        element.dispatchEvent(new Event(eventName));
        element.dispatchEvent(new MouseEvent("dblclick", { bubbles: true }));
      };
      reset(master, -30);
      reset(hardness, 0.1);
      reset(colour, 50);
      reset(mode, "difference", "change");
      reset(decay, 0.9);
      reset(model, Math.min(12, initial.model + 2));
      const deadline = performance.now() + 5000;
      const poll = () => {
        if (document.getElementById("status").textContent === "Ready" ||
            performance.now() >= deadline) {
          const checks = {
            master: master.value === "-12",
            hardness: hardness.value === initial.hardness,
            colour: colour.value === "90",
            mode: mode.value === "mirror",
            decay: decay.value === initial.decay,
            model: Math.abs(Number(model.value) - initial.model) <= 0.11,
          };
          resolve({ checks, passed: Object.values(checks).every(Boolean) });
        } else setTimeout(poll, 25);
      };
      poll();
    })`,
    awaitPromise: true, returnByValue: true,
  });
  result.controls = controls.result.value;
}
if (screenshot) {
  const capture = await call("Page.captureScreenshot", { format: "png" });
  await writeFile(screenshot, Buffer.from(capture.data, "base64"));
}
socket.close();
console.log(JSON.stringify(result));
const audioMeter = result.audio?.afterStrikes?.meter ?? "";
const xruns = Number(audioMeter.match(/xruns (\d+)/)?.[1] ?? 0);
if (result.status !== "Ready" || result.lit === 0 || testAudio &&
    (!audioMeter.includes("running") || !/hits [1-9]\d*/.test(audioMeter) ||
      !audioMeter.includes("out -") || xruns > 32) ||
    testControls && !result.controls?.passed) {
  process.exitCode = 1;
}
