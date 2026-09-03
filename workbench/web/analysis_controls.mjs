const byId = id => document.getElementById(id);

function bindDivider() {
  const divider = byId("analysis-divider");
  const analysis = document.querySelector(".analysis");
  const storageKey = "tf-spectrogram-height-v2";
  const apply = height => {
    const maximum = Math.min(620, Math.max(220, .65 * window.innerHeight));
    const value = Math.round(Math.max(110, Math.min(maximum, height)));
    analysis.style.setProperty("--spectrogram-height", `${value}px`);
    divider.setAttribute("aria-valuenow", value);
    return value;
  };
  let height = apply(Number(localStorage.getItem(storageKey)) || 300);
  let drag = null;
  divider.onpointerdown = event => {
    if (event.button !== 0) return;
    event.preventDefault();
    drag = { pointerId: event.pointerId, y: event.clientY, height };
    divider.setPointerCapture(event.pointerId);
  };
  divider.onpointermove = event => {
    if (!drag || event.pointerId !== drag.pointerId) return;
    height = apply(drag.height + event.clientY - drag.y);
  };
  const finish = event => {
    if (!drag || event.pointerId !== drag.pointerId) return;
    drag = null;
    localStorage.setItem(storageKey, height);
    if (divider.hasPointerCapture(event.pointerId))
      divider.releasePointerCapture(event.pointerId);
  };
  divider.onpointerup = finish;
  divider.onpointercancel = finish;
  divider.ondblclick = event => {
    event.preventDefault();
    height = apply(300);
    localStorage.setItem(storageKey, height);
  };
  divider.onkeydown = event => {
    if (!["ArrowUp", "ArrowDown", "Home"].includes(event.key)) return;
    event.preventDefault();
    height = apply(event.key === "Home" ? 300 :
      height + (event.key === "ArrowDown" ? 12 : -12));
    localStorage.setItem(storageKey, height);
  };
}

export function bindAnalysisControls({
  state, view, analyzeReference, analyzeSynthesis, scheduleRender,
}) {
  bindDivider();
  byId("view-mode").onchange = event =>
    view.setSettings({ mode: event.target.value });
  const update = () => {
    state.analysis.size = Number(byId("fft-size").value);
    state.analysis.hop = Math.round(state.analysis.size *
      (1 - Number(byId("overlap").value)));
    state.analysis.window = byId("window").value;
    analyzeReference();
    if (state.synthesis) analyzeSynthesis();
  };
  byId("fft-size").onchange = update;
  byId("overlap").onchange = update;
  byId("window").onchange = update;
  byId("render-seconds").onchange = scheduleRender;
  byId("colour-range").oninput = event => {
    state.analysis.dynamicRangeDb = Number(event.target.value);
    event.target.nextElementSibling.textContent = `${event.target.value} dB`;
    view.setSettings({ dynamicRangeDb: state.analysis.dynamicRangeDb });
  };
  byId("reset-view").onclick = () => view.reset();
  const resetSelect = (id, value, eventName = "change") => {
    byId(id).ondblclick = event => {
      event.preventDefault();
      event.currentTarget.value = value;
      event.currentTarget.dispatchEvent(new Event(eventName));
    };
  };
  resetSelect("view-mode", "mirror");
  resetSelect("window", "hann");
  resetSelect("fft-size", "2048");
  resetSelect("overlap", "0.75");
  resetSelect("render-seconds", "reference");
  byId("colour-range").ondblclick = event => {
    event.preventDefault();
    event.currentTarget.value = 90;
    event.currentTarget.dispatchEvent(new Event("input"));
  };
}
