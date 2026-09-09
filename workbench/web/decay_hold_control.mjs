// Design-time assistance only. Accepted edits are ordinary visible parameters.
const Keys = new Set(['bloom_rate', 'bloom_energy_acceleration',
  'body_brightness', 'body_excitation_centre', 'bloom_timing_meta']);
const signature = value => JSON.stringify(value);

export function mountDecayHold(parent, {state, read, apply, onError}) {
  const panel = document.createElement('div');
  panel.className = 'decay-hold';
  panel.innerHTML = `<label><input type="checkbox"> Hold decay</label>
    <button type="button" hidden>Cancel</button><div role="status"></div>`;
  parent.append(panel);
  panel.dataset.tooltip = 'After a Bloom edit, try to preserve the previous 1–6 second tail using the visible T60 knots and, unless you edited it, diffusion nonlinearity. The attack follows your edit. No gain matching or extra envelope. Large changes may not be compensable.';
  const toggle = panel.querySelector('input'), cancel = panel.querySelector('button');
  const status = panel.querySelector('[role=status]'), events = new AbortController();
  toggle.checked = state.holdDecayEnabled ??= true;
  let gesture, worker, timer, pending, started, evaluations = 0, expected;
  const keyOf = target => target.closest('[data-fit-key]')?.dataset.fitKey;
  function stop(message) {
    clearTimeout(pending); clearInterval(timer);
    worker?.terminate(); worker = null; gesture = null; cancel.hidden = true;
    if (message) status.textContent = message;
  }
  function start(event) {
    const key = keyOf(event.target);
    if (event.type === 'keydown' && !['ArrowLeft','ArrowRight','ArrowUp','ArrowDown','Home','End','PageUp','PageDown'].includes(event.key)) return;
    if (gesture?.target === event.target) return;
    stop(worker ? 'Cancelled by a new edit.' : '');
    if (toggle.checked && Keys.has(key) && event.target.matches('input[type=range]'))
      gesture = {target:event.target, context:read()};
  }
  function finish(event) {
    if (!gesture || gesture.target !== event.target) return;
    const baseline = gesture.context, edited = read();
    gesture = null;
    if (signature(baseline) === signature(edited)) return;
    if (signature({...baseline, parameters:[]}) !== signature({...edited, parameters:[]})) return;
    // Delay allows a double-click reset to cancel instead of fighting it.
    pending = setTimeout(() => run(baseline, edited), 350);
  }
  function run(baseline, edited) {
    if (!toggle.checked || signature(read()) !== signature(edited)) return;
    expected = signature(edited); started = performance.now(); evaluations = 0;
    cancel.hidden = false;
    const tick = () => {
      if (signature(read()) !== expected) return stop('Cancelled: sound or strike changed.');
      status.textContent = `Holding decay… ${((performance.now()-started)/1000).toFixed(1)} s · ${evaluations} renders`;
    };
    tick(); timer = setInterval(tick, 200);
    try {
      worker = new Worker('decay_hold_worker.mjs', {type:'module'});
      worker.onmessage = ({data}) => {
        if (data.progress) { evaluations = data.progress; return; }
        if (signature(read()) !== expected) return stop('Cancelled: sound changed.');
        stop();
        if (data.error) { status.textContent = data.error; onError(new Error(data.error)); return; }
        const result = data.result;
        if (result.accepted) apply(result.values);
        status.textContent = `${result.reason} · ${result.before.toFixed(2)} → ${result.after.toFixed(2)} dB tail change${result.atLimit ? ' · control limit' : ''}${result.accepted ? ' · visible controls updated' : ' · your edit kept'}`;
      };
      worker.onerror = event => { stop('Hold decay failed; your edit is unchanged.'); onError(new Error(event.message)); };
      worker.postMessage({...edited, baseline:baseline.parameters, edited:edited.parameters});
    } catch (error) { stop('Hold decay failed; your edit is unchanged.'); onError(error); }
  }
  // One anchor per drag / keyboard gesture, not one correction per input tick.
  document.addEventListener('pointerdown', start, {capture:true, signal:events.signal});
  parent.addEventListener('keydown', start, {capture:true, signal:events.signal});
  parent.addEventListener('change', finish, {signal:events.signal});
  parent.addEventListener('dblclick', () => stop('Reset applied without compensation.'), {capture:true, signal:events.signal});
  parent.addEventListener('pointercancel', () => stop('Cancelled.'), {signal:events.signal});
  toggle.onchange = () => { state.holdDecayEnabled = toggle.checked; stop(toggle.checked ? 'Ready for a Bloom edit.' : 'Off — decay controls stay fixed.'); };
  cancel.onclick = () => stop('Cancelled; your edit is unchanged.');
  status.textContent = 'Release a Bloom control to compensate · may take a few seconds';
  return {destroy() {stop(); events.abort();}};
}
