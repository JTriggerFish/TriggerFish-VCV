import { writeFile } from "node:fs/promises";

const arguments_ = process.argv.slice(2);
const endpoint = arguments_.find(argument => /^https?:\/\//.test(argument)) ??
  "http://127.0.0.1:9223";
const screenshot = arguments_.find(argument =>
  argument !== endpoint && !argument.startsWith("--"));
const testAudio = arguments_.includes("--trigger");
const profileUi = arguments_.includes("--profile");
const testControls = arguments_.includes("--controls");
const reloadPage = arguments_.includes("--reload");

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
if (reloadPage) {
  await call("Page.enable");
  await call("Page.reload", { ignoreCache: true });
  await new Promise(resolve => setTimeout(resolve, 500));
}
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
      const implementChoices = [...document.querySelectorAll(
        'input[name="implement"]',
      )];
      const sizeMeta = document.getElementById("size-meta");
      const colour = document.getElementById("colour-range");
      const mode = document.getElementById("view-mode");
      const model = document.querySelector("#model-level input");
      const shape = document.querySelector("#body-controls input");
      const washDensity = document.querySelector(
        '[data-fit-key="dense_mode_density"] input[type="range"]');
      const bodyLowT60 = document.querySelector(
        '[data-fit-key="body_decay_seconds_0"] input');
      const sidebar = document.querySelector("aside");
      const analysis = document.querySelector(".analysis");
      const analysisTop = analysis.getBoundingClientRect().top;
      const spectrogram = document.getElementById("spectrogram");
      const referencePixels = () => spectrogram.getContext("2d").getImageData(
        0, 0, Math.floor(spectrogram.width / 2), spectrogram.height,
      ).data;
      const referenceBefore = referencePixels();
      sidebar.scrollTop = Math.min(240, sidebar.scrollHeight - sidebar.clientHeight);
      const initial = {
        hardness: hardness.value,
        implement: implementChoices.find(input => input.checked)?.value,
        model: Number(model.value), shape: shape.value,
      };
      const reset = (element, changed, eventName = "input") => {
        element.value = changed;
        element.dispatchEvent(new Event(eventName));
        element.dispatchEvent(new MouseEvent("dblclick", { bubbles: true }));
      };
      reset(master, -30);
      reset(hardness, 0.1);
      implementChoices.find(input => input.value === "0").click();
      const brushCharacter = document.getElementById("character-label").textContent;
      implementChoices.find(input => input.value === "0.5").click();
      const malletCharacter = document.getElementById("character-label").textContent;
      implementChoices.find(input => input.value === initial.implement).click();
      reset(sizeMeta, 0.1);
      reset(colour, 50);
      reset(mode, "difference", "change");
      reset(shape, 0.9);
      reset(model, Math.min(1, initial.model + 0.1));
      const deadline = performance.now() + 5000;
      const poll = () => {
        if (document.getElementById("status").textContent === "Ready" ||
            performance.now() >= deadline) {
          const checks = {
            master: master.value === "-12",
            hardness: hardness.value === initial.hardness,
            implement: implementChoices.find(input => input.checked)?.value ===
              initial.implement,
            implementChoices: implementChoices.length === 3,
            contextualCharacter: document.getElementById("character-label")
              .textContent === ({
                "0": "Bristle stiffness", "0.5": "Mallet firmness",
                "1": "Tip hardness",
              })[initial.implement],
            familyCharacterLabels: brushCharacter === "Bristle stiffness" &&
              malletCharacter === "Mallet firmness",
            sizeMeta: sizeMeta.value === "0.5",
            colour: colour.value === "90",
            mode: mode.value === "mirror",
            shape: shape.value === initial.shape,
            model: Math.abs(Number(model.value) - initial.model) <= 0.11,
            bodyLowT60: bodyLowT60 instanceof HTMLInputElement,
            fiveBodyT60Knots: document.querySelectorAll(
              "#decay-editor .editor-point").length === 5,
            sharedT60Editor: Boolean(document.querySelector("#decay-editor svg")),
            resolvedModeBars: document.querySelectorAll(
              "#resolved-editor .editor-bar")
              .length === 12,
            resolvedPaintMode: Boolean(document.querySelector(
              "#resolved-editor .editor-mode input[type=checkbox]",
            )),
            noResolvedModeToggle: !document.querySelector(
              '[data-fit-key="resolved_modes_enabled"]'),
            denseWashCurve: document.querySelectorAll(
              "#dense-wash-editor .editor-point").length === 8,
            continuousWashDensity: washDensity instanceof HTMLInputElement &&
              washDensity.step !== "1",
            waveformDragDisabled:
              document.getElementById("waveform")._fullLayout?.dragmode === false,
            turbulenceDefaultOn: document.querySelector(
              "#turbulence-toggle input")?.checked === true,
            noStrikeEllipse: getComputedStyle(
              document.getElementById("strike-pad"), "::after").content === "none",
            twoControlColumns: getComputedStyle(
              document.querySelector(".control-columns"),
            ).gridTemplateColumns.split(" ").length === 2,
            independentScroll: sidebar.scrollHeight > sidebar.clientHeight &&
              analysis.getBoundingClientRect().top === analysisTop,
            strikeInAnalysis: Boolean(
              document.getElementById("strike-pad").closest(".analysis")),
            strikeVisible: document.getElementById("strike-pad")
              .getBoundingClientRect().bottom <=
              analysis.getBoundingClientRect().bottom,
            turbulenceInLeftColumn: Boolean(document.querySelector(
              ".control-column:first-child #turbulence-toggle",
            )),
            referenceColourScale: document.getElementById("colour-ceiling")
              .textContent.startsWith("Reference ceiling"),
            livePreparationVisible: document.getElementById("live-commit")
              .textContent.startsWith("Live DSP"),
            fixedReferenceColours: (() => {
              const after = referencePixels();
              return after.length === referenceBefore.length && after.every(
                (value, index) => value === referenceBefore[index],
              );
            })(),
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
