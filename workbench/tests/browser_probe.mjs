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
const captureStart = process.env.TF_WORKBENCH_CAPTURE_START;

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
if (result.status === "Loading DSP…") {
  const diagnostic = await call("Runtime.evaluate", {
    expression: `import("./app.mjs?diagnostic=1")
      .then(() => ({ loaded: true }))
      .catch(error => ({ loaded: false, message: String(error), stack: error.stack }))`,
    awaitPromise: true, returnByValue: true,
  });
  result.moduleDiagnostic = diagnostic.result.value;
}
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
  await call("Runtime.evaluate", {
    expression: `new Promise(resolve => {
      const deadline = performance.now() + 15000;
      let readySince = 0;
      const poll = () => {
        const ready = document.getElementById("status")?.textContent === "Ready" &&
          document.getElementById("live-commit")?.textContent
            .startsWith("Live DSP ready");
        readySince = ready ? (readySince || performance.now()) : 0;
        if (readySince && performance.now() - readySince >= 500 ||
            performance.now() >= deadline) {
          document.getElementById("stop").click();
          resolve();
        } else setTimeout(poll, 50);
      };
      poll();
    })`,
    awaitPromise: true,
  });
  await new Promise(resolve => setTimeout(resolve, 100));
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
      const recipe = document.getElementById("instrument-recipe");
      if (recipe.value !== "0") {
        recipe.value = "0";
        recipe.dispatchEvent(new Event("change"));
      }
      const master = document.getElementById("master");
      const hardness = document.getElementById("hardness");
      const contactSpread = document.getElementById("contact-spread");
      const implementChoices = [...document.querySelectorAll(
        'input[name="implement"]',
      )];
      const sizeMeta = document.getElementById("size-meta");
      const colour = document.getElementById("colour-range");
      const mode = document.getElementById("view-mode");
      const model = document.querySelector("#model-level input");
      const bodyLowT60 = document.querySelector(
        '[data-fit-key="body_decay_seconds_0"] input');
      const fieldTurbulence = document.querySelector(
        '[data-fit-key="field_turbulence"] input');
      const settingsDialog = document.getElementById("settings-dialog");
      const routingCompact = document.getElementById("routing-compact");
      const routingDialog = document.getElementById("routing-dialog");
      routingCompact.dispatchEvent(new MouseEvent("dblclick", { bubbles: true }));
      const routingOpened = routingDialog.open;
      const routingNodeCount = document.querySelectorAll(
        "#routing-expanded .routing-node").length;
      const routingEdgeCount = document.querySelectorAll(
        "#routing-expanded .routing-edge").length;
      document.querySelector(
        '#routing-expanded [data-connection-id="route-1"]',
      ).dispatchEvent(new MouseEvent("click", { bubbles: true }));
      const routingDisabled = document.querySelector(
        '#routing-expanded [data-connection-id="route-1"]',
      ).classList.contains("disabled");
      document.querySelector(
        '#routing-expanded [data-connection-id="route-1"]',
      ).dispatchEvent(new MouseEvent("click", { bubbles: true }));
      const routingRestored = !document.querySelector(
        '#routing-expanded [data-connection-id="route-1"]',
      ).classList.contains("disabled");
      document.getElementById("routing-close").click();
      const routingClosed = !routingDialog.open;
      document.getElementById("settings-open").click();
      const settingsOpened = settingsDialog.open;
      const directWorklet = document.getElementById("audio-buffer-status")
        .textContent.includes("128-frame direct worklet");
      const midiSettingsPresent =
        document.getElementById("midi-enable") instanceof HTMLButtonElement &&
        document.getElementById("midi-input") instanceof HTMLSelectElement &&
        document.getElementById("midi-channel") instanceof HTMLSelectElement;
      const midiEnabledByDefault =
        document.getElementById("midi-status").textContent !==
          "MIDI access has not been requested.";
      document.getElementById("settings-close").click();
      const settingsClosed = !settingsDialog.open;
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
        contactSpread: contactSpread.value,
        implement: implementChoices.find(input => input.checked)?.value,
        model: Number(model.value),
      };
      const reset = (element, changed, eventName = "input") => {
        element.value = changed;
        element.dispatchEvent(new Event(eventName));
        element.dispatchEvent(new MouseEvent("dblclick", { bubbles: true }));
      };
      reset(master, -30);
      reset(hardness, 0.1);
      reset(contactSpread, 0.9);
      implementChoices.find(input => input.value === "0").click();
      const brushCharacter = document.getElementById("character-label").textContent;
      implementChoices.find(input => input.value === "0.5").click();
      const malletCharacter = document.getElementById("character-label").textContent;
      implementChoices.find(input => input.value === initial.implement).click();
      reset(sizeMeta, 0.1);
      reset(colour, 50);
      reset(mode, "difference", "change");
      reset(model, Math.min(1, initial.model + 0.1));
      const modal = document.querySelector("#modal-editor svg");
      const modalCount = modal.querySelectorAll(".modal-bar").length;
      modal.querySelector(".modal-bar").dispatchEvent(new MouseEvent(
        "dblclick", { bubbles: true },
      ));
      const modalDeleteWorked =
        modal.querySelectorAll(".modal-bar").length === modalCount - 1;
      const modalPreset = document.getElementById("modal-preset");
      modalPreset.value = "fitted";
      modalPreset.dispatchEvent(new Event("change"));
      const modalRestoreWorked =
        modal.querySelectorAll(".modal-bar").length === modalCount;
      const modalSelection = document.getElementById("modal-selection");
      const persistentModalInspector =
        modalSelection.querySelectorAll('input[type="range"]').length === 3 &&
        [...modalSelection.querySelectorAll('input[type="range"]')]
          .every(input => input.disabled) &&
        modalSelection.getClientRects().length > 0;
      const harmonicGuide = document.getElementById("harmonic-guide");
      harmonicGuide.checked = true;
      harmonicGuide.dispatchEvent(new Event("change"));
      const harmonicGuideWorks =
        modal.querySelectorAll(".harmonic-grid").length > 1 &&
        !document.getElementById("harmonic-note").disabled &&
        !document.getElementById("harmonic-snap").disabled &&
        !document.getElementById("harmonic-snap-all").disabled;
      const decayEditor = document.querySelector("#decay-editor svg");
      const decayCount = decayEditor.querySelectorAll(".editor-point").length;
      const boundaryDelete = document.querySelector("#decay-selection button");
      const boundaryDeleteIsClear = boundaryDelete?.textContent === "Delete knot" &&
        boundaryDelete.disabled;
      const decayBounds = decayEditor.getBoundingClientRect();
      decayEditor.dispatchEvent(new MouseEvent("dblclick", {
        bubbles: true,
        clientX: decayBounds.left + .63 * decayBounds.width,
        clientY: decayBounds.top + .45 * decayBounds.height,
      }));
      const decayInsertWorked =
        decayEditor.querySelectorAll(".editor-point").length === decayCount + 1;
      const interiorDelete = document.querySelector("#decay-selection button");
      const interiorDeleteIsEnabled = interiorDelete?.textContent === "Delete knot" &&
        !interiorDelete.disabled;
      interiorDelete?.click();
      const decayRemoveWorked =
        decayEditor.querySelectorAll(".editor-point").length === decayCount;
      const divider = document.getElementById("analysis-divider");
      const dividerStart = divider.getAttribute("aria-valuenow");
      divider.dispatchEvent(new KeyboardEvent("keydown", {
        key: "ArrowDown", bubbles: true,
      }));
      divider.dispatchEvent(new KeyboardEvent("keydown", {
        key: "ArrowUp", bubbles: true,
      }));
      const modalDesign = document.querySelector(".modal-design");
      const modalBounds = modal.getBoundingClientRect();
      const modalDesignBounds = modalDesign.getBoundingClientRect();
      const modalAspect = modal.viewBox.baseVal.width / modal.viewBox.baseVal.height;
      const deadline = performance.now() + 5000;
      const poll = () => {
        if (document.getElementById("status").textContent === "Ready" ||
            performance.now() >= deadline) {
          const checks = {
            modularRouting: routingOpened && routingClosed &&
              routingNodeCount === 5 && routingEdgeCount === 6 &&
              routingDisabled && routingRestored &&
              document.querySelectorAll("[data-module-id]").length >= 8 &&
              Boolean(document.querySelector(
                '[data-module-id="body"][style*="--module-colour"]')),
            master: master.value === "-12",
            calibrationPicker: (() => {
              const picker = document.getElementById("instrument-calibration");
              const labels = [...picker.options].map(option => option.textContent);
              return picker.closest(".calibration-picker") &&
                labels.some(label => label.startsWith("Snare — medium standard hit")) &&
                labels.some(label => label.startsWith("Acoustic kick — medium centre")) &&
                labels.some(label => label.startsWith("Gong — representative mallet")) &&
                labels.some(label => label.startsWith("Ride — medium bow"));
            })(),
            settingsMenu: settingsOpened && settingsClosed && directWorklet &&
              midiSettingsPresent && midiEnabledByDefault,
            hardness: hardness.value === initial.hardness,
            contactSpread: contactSpread.value === initial.contactSpread,
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
            model: Math.abs(Number(model.value) - initial.model) <= 0.11,
            bodyLowT60: bodyLowT60 instanceof HTMLInputElement,
            flexibleBodyT60Knots:
              document.querySelectorAll("#decay-editor .editor-point").length === 2 &&
              document.querySelectorAll("#decay-editor rect.editor-point").length === 2 &&
              Boolean(document.querySelector("#decay-editor .decay-all-handle")) &&
              boundaryDeleteIsClear && interiorDeleteIsEnabled &&
              decayInsertWorked && decayRemoveWorked,
            sharedT60Editor: Boolean(decayEditor),
            unifiedTurbulenceControl:
              fieldTurbulence instanceof HTMLInputElement,
            noLegacyBodyUi: !document.getElementById("body-ui-mode") &&
              !document.querySelector('[data-ui-mode="legacy"]'),
            bloomDiffusionControl: Boolean(document.querySelector(
              '[data-fit-key="bloom_diffusion"] input')),
            independentBloomControls: Boolean(document.querySelector(
              '[data-fit-key="bloom_level"] input')) &&
              Boolean(document.querySelector(
                '[data-fit-key="bloom_nonlinearity"] input')),
            modalPacketEditor: Boolean(document.querySelector(
              "#modal-editor svg.modal-editor")) &&
              document.querySelectorAll("#modal-editor .modal-bar").length > 0 &&
              document.querySelectorAll("#modal-editor .modal-packet").length > 0,
            modalEditingTools: Boolean(document.getElementById("modal-tool-edit")) &&
              Boolean(document.getElementById("modal-tool-level")) &&
              Boolean(document.getElementById("modal-tool-paint")) &&
              Boolean(document.getElementById("modal-clear")) &&
              modalDeleteWorked && modalRestoreWorked,
            persistentModalInspector,
            harmonicGuideWorks,
            modalCeiling15k: [...document.querySelectorAll(
              "#modal-editor .editor-tick")].some(node =>
                node.textContent === "15k") &&
              ![...document.querySelectorAll("#modal-editor .editor-tick")]
                .some(node => node.textContent === "20k"),
            modalTooltips: Boolean(document.querySelector(
              "#modal-tool-edit[data-tooltip]")) && Boolean(document.querySelector(
                '[data-fit-key="field_turbulence"][data-tooltip]')),
            modalUsesAvailableSpace:
              modalBounds.width >= modalDesignBounds.width - 30 &&
              Math.abs(modalAspect - modalBounds.width / modalBounds.height) < .05 &&
              modalDesignBounds.height >= 375 &&
              modalBounds.bottom <= modalDesignBounds.bottom &&
              spectrogram.getBoundingClientRect().bottom <
                document.getElementById("strike-pad").getBoundingClientRect().top &&
              document.getElementById("strike-pad").getBoundingClientRect().bottom <
                modalDesignBounds.top,
            noResolvedModeToggle: !document.querySelector(
              '[data-fit-key="resolved_modes_enabled"]'),
            waveformDragDisabled:
              document.getElementById("waveform")._fullLayout?.dragmode === false,
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
            resizableSpectrogram: document.getElementById("analysis-divider")
              .getAttribute("role") === "separator" &&
              divider.getAttribute("aria-valuenow") === dividerStart,
            compactRadiation: Boolean(document.getElementById(
              "direct-radiation-advanced")) && Boolean(document.getElementById(
                "dense-radiation-advanced")),
          };
          resolve({
            checks, passed: Object.values(checks).every(Boolean),
            liveStatus: document.getElementById("live-commit").textContent,
            layout: {
              modalWidth: modalBounds.width,
              designWidth: modalDesignBounds.width,
              modalBottom: modalBounds.bottom,
              designBottom: modalDesignBounds.bottom,
              analysisBottom: analysis.getBoundingClientRect().bottom,
              viewBoxAspect: modalAspect,
              elementAspect: modalBounds.width / modalBounds.height,
            },
          });
        } else setTimeout(poll, 25);
      };
      poll();
    })`,
    awaitPromise: true, returnByValue: true,
  });
  if (!controls.result.value?.checks) {
    throw new Error(`metallic control probe failed: ${
      controls.exceptionDetails?.exception?.description ??
      controls.exceptionDetails?.text ?? JSON.stringify(controls)}`);
  }
  result.controls = controls.result.value;
  const kickRecipe = await call("Runtime.evaluate", {
    expression: `new Promise(resolve => {
      const select = document.getElementById("instrument-recipe");
      select.value = "1";
      select.dispatchEvent(new Event("change"));
      const deadline = performance.now() + 12000;
      const poll = () => {
        const ready = document.getElementById("status").textContent === "Ready" &&
          document.getElementById("live-commit").textContent
            .startsWith("Live DSP ready");
        if (!ready && performance.now() < deadline) {
          setTimeout(poll, 50);
          return;
        }
        document.getElementById("routing-compact").dispatchEvent(
          new MouseEvent("dblclick", { bubbles: true }),
        );
        const firstRoute = document.querySelector(
          '#routing-expanded [data-connection-id="kick-route-1"]',
        );
        firstRoute?.dispatchEvent(new MouseEvent("click", { bubbles: true }));
        const routeDisabled = document.querySelector(
          '#routing-expanded [data-connection-id="kick-route-1"]',
        )?.classList.contains("disabled");
        document.querySelector(
          '#routing-expanded [data-connection-id="kick-route-1"]',
        )?.dispatchEvent(new MouseEvent("click", { bubbles: true }));
        document.getElementById("routing-close").click();
        const hitsBefore = Number(document.getElementById("limiter").textContent
          .match(/hits (\\d+)/)?.[1] ?? 0);
        document.getElementById("play-synthesis").click();
        let loudestDb = -Infinity;
        const sampleMeter = setInterval(() => {
          const match = document.getElementById("limiter").textContent
            .match(/out (-?\\d+) dBFS/);
          if (match) loudestDb = Math.max(loudestDb, Number(match[1]));
        }, 25);
        setTimeout(() => {
          clearInterval(sampleMeter);
          const hitsAfter = Number(document.getElementById("limiter").textContent
            .match(/hits (\\d+)/)?.[1] ?? 0);
          const checks = {
            selector: select.options.length === 4 && select.value === "1",
            controls: document.querySelectorAll(
              '[data-kick-key] input[type="range"]',
            ).length === 15,
            contextualPanels:
              getComputedStyle(document.querySelector(
                '[data-module-id="kick-primary"]')).display !== "none" &&
              [...document.querySelectorAll('[data-module-id="body"]')]
                .every(element => getComputedStyle(element).display === "none"),
            routing: document.querySelectorAll(
              "#routing-compact .routing-node",
            ).length === 6 && document.querySelectorAll(
              "#routing-compact .routing-edge",
            ).length === 5 && routeDisabled,
            label: document.getElementById("routing-recipe-label")
              .textContent === "Compact FM kick",
            liveAudio: hitsAfter > hitsBefore && loudestDb > -80,
          };
          resolve({
            checks, passed: Object.values(checks).every(Boolean),
            liveStatus: document.getElementById("live-commit").textContent,
            loudestDb,
          });
        }, 500);
      };
      poll();
    })`,
    awaitPromise: true, returnByValue: true,
  });
  result.controls.kickRecipe = kickRecipe.result.value;
  result.controls.checks.kickRecipe = result.controls.kickRecipe.passed;
  const membraneRecipe = await call("Runtime.evaluate", {
    expression: `new Promise(resolve => {
      const select = document.getElementById("instrument-recipe");
      select.value = "2";
      select.dispatchEvent(new Event("change"));
      const deadline = performance.now() + 12000;
      const poll = () => {
        const ready = document.getElementById("status").textContent === "Ready" &&
          document.getElementById("live-commit").textContent
            .startsWith("Live DSP ready");
        if (!ready && performance.now() < deadline) {
          setTimeout(poll, 50);
          return;
        }
        const preset = [...document.querySelectorAll(
          "#membrane-preset-controls button",
        )].find(button => button.textContent === "Acoustic kick");
        preset?.click();
        const hitsBefore = Number(document.getElementById("limiter").textContent
          .match(/hits (\\d+)/)?.[1] ?? 0);
        document.getElementById("play-synthesis").click();
        let loudestDb = -Infinity;
        const sampleMeter = setInterval(() => {
          const match = document.getElementById("limiter").textContent
            .match(/out (-?\\d+) dBFS/);
          if (match) loudestDb = Math.max(loudestDb, Number(match[1]));
        }, 25);
        setTimeout(() => {
          clearInterval(sampleMeter);
          const hitsAfter = Number(document.getElementById("limiter").textContent
            .match(/hits (\\d+)/)?.[1] ?? 0);
          const checks = {
            selector: select.options.length === 4 && select.value === "2",
            controls: document.querySelectorAll(
              '[data-membrane-key] input[type="range"]',
            ).length === 32,
            presets: Boolean(preset) && document.querySelector(
              '[data-membrane-key="fundamental_hz"] output',
            )?.textContent.startsWith("35.0"),
            equalizerModes: document.querySelectorAll(
              "#membrane-eq-mode-controls option",
            ).length === 3,
            contextualPanels: getComputedStyle(document.querySelector(
              '[data-module-id="membrane-body"]',
            )).display !== "none" && getComputedStyle(document.querySelector(
              '[data-module-id="kick-primary"]',
            )).display === "none",
            routing: document.querySelectorAll(
              "#routing-compact .routing-node",
            ).length === 9 && document.querySelectorAll(
              "#routing-compact .routing-edge",
            ).length === 11,
            label: document.getElementById("routing-recipe-label")
              .textContent === "Membrane drum",
            liveAudio: hitsAfter > hitsBefore && loudestDb > -80,
          };
          resolve({ checks, passed: Object.values(checks).every(Boolean), loudestDb });
        }, 500);
      };
      poll();
    })`,
    awaitPromise: true, returnByValue: true,
  });
  result.controls.membraneRecipe = membraneRecipe.result.value;
  result.controls.checks.membraneRecipe =
    result.controls.membraneRecipe.passed;
  const snareRecipe = await call("Runtime.evaluate", {
    expression: `new Promise(resolve => {
      const select = document.getElementById("instrument-recipe");
      select.value = "3";
      select.dispatchEvent(new Event("change"));
      const deadline = performance.now() + 12000;
      const poll = () => {
        const ready = document.getElementById("status").textContent === "Ready" &&
          document.getElementById("live-commit").textContent
            .startsWith("Live DSP ready");
        if (!ready && performance.now() < deadline) {
          setTimeout(poll, 50);
          return;
        }
        const hitsBefore = Number(document.getElementById("limiter").textContent
          .match(/hits (\\d+)/)?.[1] ?? 0);
        document.getElementById("play-synthesis").click();
        setTimeout(() => {
          const hitsAfter = Number(document.getElementById("limiter").textContent
            .match(/hits (\\d+)/)?.[1] ?? 0);
          const checks = {
            selector: select.options.length === 4 && select.value === "3",
            wireControls: document.querySelectorAll(
              "#snare-wire-response-controls input, #snare-wire-spectrum-controls input",
            ).length === 14,
            ringControls: document.querySelectorAll(
              "#snare-ring-controls input",
            ).length === 3,
            routing: document.querySelectorAll(
              "#routing-compact .routing-node",
            ).length === 10 && document.querySelectorAll(
              "#routing-compact .routing-edge",
            ).length === 13,
            label: document.getElementById("routing-recipe-label")
              .textContent === "Snare drum",
            liveAudio: hitsAfter > hitsBefore,
          };
          resolve({ checks, passed: Object.values(checks).every(Boolean) });
        }, 500);
      };
      poll();
    })`,
    awaitPromise: true, returnByValue: true,
  });
  result.controls.snareRecipe = snareRecipe.result.value;
  result.controls.checks.snareRecipe = result.controls.snareRecipe.passed;
  const calibrationPreset = await call("Runtime.evaluate", {
    expression: `new Promise(resolve => {
      const picker = document.getElementById("instrument-calibration");
      picker.value = "snare-standard";
      picker.dispatchEvent(new Event("change"));
      const deadline = performance.now() + 12000;
      const poll = () => {
        const ready = document.getElementById("status").textContent === "Ready";
        const loaded = picker.value === "snare-standard" &&
          document.getElementById("instrument-recipe").value === "3" &&
          document.getElementById("reference-corpus").value ===
            "acoustic-snare-maple" &&
          document.getElementById("reference-articulation").value === "main" &&
          document.getElementById("reference-velocity").nextElementSibling
            .textContent === "v082";
        if (ready && loaded || performance.now() >= deadline) {
          resolve({
            passed: ready && loaded &&
              document.querySelector(
                '[data-membrane-key="fundamental_hz"] output',
              )?.textContent.startsWith("185.0") &&
              document.querySelector(
                '[data-membrane-key="model_level_db"] output',
              )?.textContent.startsWith("-10.00") &&
              document.querySelector(
                'input[name="membrane-implement"][value="1"]',
              )?.checked,
            ready, loaded,
          });
        } else setTimeout(poll, 50);
      };
      poll();
    })`,
    awaitPromise: true, returnByValue: true,
  });
  result.controls.calibrationPreset = calibrationPreset.result.value;
  result.controls.checks.calibrationPreset =
    result.controls.calibrationPreset.passed;
  result.controls.passed = Object.values(result.controls.checks).every(Boolean);
}
if (screenshot) {
  if (captureStart) {
    const selected = await call("Runtime.evaluate", {
      expression: `new Promise(resolve => {
        const picker = document.getElementById("instrument-calibration");
        picker.value = ${JSON.stringify(captureStart)};
        picker.dispatchEvent(new Event("change"));
        const deadline = performance.now() + 15000;
        const poll = () => {
          const ready = document.getElementById("status").textContent === "Ready";
          if (ready || performance.now() >= deadline) resolve({ ready });
          else setTimeout(poll, 50);
        };
        poll();
      })`,
      awaitPromise: true, returnByValue: true,
    });
    result.captureStart = { id: captureStart, ...selected.result.value };
  }
  await call("Runtime.evaluate", {
    expression: `(() => {
      document.querySelector("aside").scrollTop = 0;
      document.querySelector(".analysis").scrollTop = 0;
    })()`,
  });
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
