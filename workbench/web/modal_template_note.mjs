// A design-time shortcut for the explicit Hz field, not a hidden pitch rule.
export function noteFrequency(note, octave) {
  return 440 * 2 ** ((12 * (octave + 1) + note - 69) / 12);
}

export function mountTemplateNote(parent, frequency) {
  const row = document.createElement("label");
  row.textContent = "Base note ";
  const note = document.createElement("select");
  note.setAttribute("aria-label", "Generator base note");
  note.append(new Option("Custom Hz", "custom"));
  ["C", "C♯", "D", "D♯", "E", "F", "F♯", "G", "G♯", "A", "A♯", "B"]
    .forEach((name, index) => note.append(new Option(name, index)));
  const octave = document.createElement("select");
  octave.setAttribute("aria-label", "Generator base octave");
  for (let value = 0; value <= 9; ++value) octave.append(new Option(value, value));
  const sync = () => {
    const hz = Number(frequency.value);
    const midi = Math.round(69 + 12 * Math.log2(hz / 440));
    const oct = Math.floor(midi / 12) - 1;
    const pitch = ((midi % 12) + 12) % 12;
    const exact = hz > 0 && oct >= 0 && oct <= 9 &&
      Math.abs(hz - noteFrequency(pitch, oct)) < .001;
    note.value = exact ? String(pitch) : "custom";
    if (exact) octave.value = String(oct);
    octave.disabled = !exact;
  };
  const choose = () => {
    octave.disabled = note.value === "custom";
    if (!octave.disabled) {
      frequency.value = noteFrequency(Number(note.value), Number(octave.value));
      frequency.dispatchEvent(new Event("input", {bubbles:true}));
    }
  };
  note.onchange = choose; octave.onchange = choose;
  frequency.addEventListener("input", sync);
  row.append(note, octave); parent.insertBefore(row, frequency.parentElement); sync();
  return sync;
}
