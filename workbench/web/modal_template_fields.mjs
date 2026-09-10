// Compact paired sliders / numeric entry for the design-time generator.
export function templateField(parent, key, title, initial, min, max, step, slider = true) {
  const row = document.createElement("label");
  row.className = "template-field";
  const label = document.createElement("span"); label.textContent = title;
  const input = document.createElement("input");
  Object.assign(input, {type:"number", value:initial, min, max, step});
  input.required = true; input.dataset.templateKey = key;
  input.setAttribute("aria-label", title);
  row.append(label);
  if (slider) {
    const range = document.createElement("input");
    Object.assign(range, {type:"range", value:initial, min, max, step});
    range.setAttribute("aria-label", title + " slider");
    range.oninput = () => {
      input.value = range.value; input.dispatchEvent(new Event("input", {bubbles:true}));
    };
    input.addEventListener("input", () => { range.value = input.value; });
    row.append(range);
    input.setMaximum = value => {
      input.max = value; range.max = Math.max(min, value);
      range.value = input.value;
    };
  }
  row.append(input); parent.append(row);
  const reset = () => {
    input.value = Math.min(Number(input.max), Math.max(Number(input.min), initial));
    input.dispatchEvent(new Event("input", {bubbles:true}));
  };
  row.ondblclick = reset;
  input.ondblclick = event => { event.stopPropagation(); reset(); };
  return input;
}
