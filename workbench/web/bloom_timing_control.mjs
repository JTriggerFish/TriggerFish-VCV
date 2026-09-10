import {BloomTimingKeys, bloomTimingValues} from "./bloom_timing_meta.mjs";
import {createTimingView, placeTimingPopover} from "./bloom_timing_view.mjs";

// Captured centre and popup visibility are UI state, never DSP parameters.
export function mountBloomTiming(parent, {read, apply, descriptors}) {
  const {button, panel, input, output, status, paint} = createTimingView(parent, descriptors);
  const events = new AbortController();
  const current = () => Object.fromEntries(BloomTimingKeys.map(key => [key, read(key)]));
  let baseline, frame;
  const highlight = active => BloomTimingKeys.forEach(key =>
    document.querySelector('[data-fit-key="'+key+'"]')?.classList.toggle("bloom-meta-target", active));
  const rebase = () => {
    baseline = current(); input.value = 0; output.textContent = "Centre";
    input.disabled = baseline.bloom_rate === 0;
    panel.querySelector('[data-action="reset"]').disabled = input.disabled;
    status.textContent = input.disabled ? "Enable Diffusion strength to use timing." : "";
    paint(baseline, baseline);
  };
  const change = () => {
    const position = Number(input.value);
    const {values, limited} = bloomTimingValues(baseline, position, descriptors);
    apply(values); paint(baseline, current());
    output.textContent = position === 0 ? "Centre" :
      (position > 0 ? "Later" : "Earlier")+" "+Math.abs(position).toFixed(2);
    status.textContent = limited.length ? "At control limit: "+limited.join(", ") : "";
    if (panel.matches(":popover-open")) placeTimingPopover(button, panel);
  };
  input.oninput = change;
  input.ondblclick = event => { event.preventDefault(); input.value = 0; change(); };
  panel.querySelector('[data-action="centre"]').onclick = rebase;
  panel.querySelector('[data-action="centre"]').title = "Use the current sound as the middle; no sound change.";
  panel.querySelector('[data-action="reset"]').onclick = () =>
    input.dispatchEvent(new MouseEvent("dblclick", {bubbles:true}));
  panel.querySelector('[data-action="reset"]').title = "Restore these three controls. Any previous Hold decay adjustments are kept.";
  panel.querySelector('[data-action="close"]').onclick = () => { panel.hidePopover(); button.focus(); };
  panel.addEventListener("beforetoggle", event => {
    const open = event.newState === "open";
    button.setAttribute("aria-expanded", String(open)); highlight(open);
    if (open) {
      paint(baseline, current());
      frame = requestAnimationFrame(() => {
        if (!panel.isConnected || !panel.matches(":popover-open")) return;
        placeTimingPopover(button, panel); input.focus({preventScroll:true});
      });
    }
  }, {signal:events.signal});
  const reposition = () => {
    if (panel.matches(":popover-open")) placeTimingPopover(button, panel);
  };
  window.addEventListener("resize", reposition, {signal:events.signal});
  document.addEventListener("scroll", reposition, {capture:true, signal:events.signal});
  rebase();
  return {rebase, destroy() {
    cancelAnimationFrame(frame); events.abort(); highlight(false);
    panel.remove(); button.remove();
  }};
}
