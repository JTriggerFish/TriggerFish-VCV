const Svg = "http://www.w3.org/2000/svg";
const View = { width: 600, height: 220, left: 48, right: 18, top: 18, bottom: 34 };
const clamp = (value, minimum, maximum) =>
  Math.max(minimum, Math.min(maximum, value));
const erb = frequency => 21.4 * Math.log10(1 + .00437 * frequency);
const inverseErb = rate => (10 ** (rate / 21.4) - 1) / .00437;

function element(name, attributes = {}) {
  const result = document.createElementNS(Svg, name);
  for (const [key, value] of Object.entries(attributes)) {
    result.setAttribute(key, value);
  }
  return result;
}

export class DecayCurveEditor {
  constructor(parent, options) {
    this.options = options;
    this.selected = 0;
    this.drag = null;
    this.svg = element("svg", {
      class: "decay-curve-editor", viewBox: `0 0 ${View.width} ${View.height}`,
      role: "application", tabindex: 0,
      "aria-label": "Frequency-dependent body T60 editor",
    });
    parent.append(this.svg);
    this.bind();
    this.paint();
  }

  nyquist() {
    return Math.max(4000, Number(this.options.nyquist()));
  }

  xPosition(frequency) {
    const amount = erb(clamp(frequency, 0, this.nyquist())) /
      erb(this.nyquist());
    return View.left + amount * (View.width - View.left - View.right);
  }

  frequency(position) {
    const amount = clamp(
      (position - View.left) / (View.width - View.left - View.right), 0, 1,
    );
    return inverseErb(amount * erb(this.nyquist()));
  }

  yPosition(logSeconds) {
    const amount = (logSeconds - this.options.minimumLogSeconds) /
      (this.options.maximumLogSeconds - this.options.minimumLogSeconds);
    return View.height - View.bottom - amount *
      (View.height - View.top - View.bottom);
  }

  logSeconds(position) {
    const amount = clamp(
      (View.height - View.bottom - position) /
        (View.height - View.top - View.bottom), 0, 1,
    );
    return this.options.minimumLogSeconds + amount *
      (this.options.maximumLogSeconds - this.options.minimumLogSeconds);
  }

  eventPosition(event) {
    const bounds = this.svg.getBoundingClientRect();
    return {
      x: View.width * (event.clientX - bounds.left) / bounds.width,
      y: View.height * (event.clientY - bounds.top) / bounds.height,
    };
  }

  constrainFrequency(slot, frequency) {
    const points = this.options.points();
    const index = points.findIndex(point => point.slot === slot);
    if (index < 0 || points[index].fixed) return points[index]?.x ?? frequency;
    const gap = .22;
    const first = erb(points[index - 1].x) + gap;
    const last = erb(points[index + 1].x) - gap;
    return inverseErb(clamp(erb(frequency), first, Math.max(first, last)));
  }

  bind() {
    this.svg.addEventListener("pointermove", event => {
      if (!this.drag || event.pointerId !== this.drag.pointerId) return;
      const position = this.eventPosition(event);
      if (this.drag.kind === "all") this.dragAll(position);
      else this.dragPoint(position);
    });
    const finish = event => {
      if (!this.drag || event.pointerId !== this.drag.pointerId) return;
      if (this.svg.hasPointerCapture(event.pointerId)) {
        this.svg.releasePointerCapture(event.pointerId);
      }
      this.drag = null;
      this.options.select(this.selected);
    };
    this.svg.addEventListener("pointerup", finish);
    this.svg.addEventListener("pointercancel", finish);
    this.svg.addEventListener("dblclick", event => {
      if (event.target.closest(".decay-handle")) return;
      event.preventDefault();
      const position = this.eventPosition(event);
      const slot = this.options.insert(
        this.frequency(position.x), this.logSeconds(position.y),
      );
      if (slot !== null) this.select(slot);
      this.paint();
    });
    this.svg.addEventListener("keydown", event => {
      if ((event.key === "Delete" || event.key === "Backspace") &&
          this.selected !== null) {
        event.preventDefault();
        this.options.remove(this.selected);
        this.selected = 0;
        this.options.select(this.selected);
        this.paint();
      }
    });
  }

  beginPoint(event, point) {
    if (event.button !== 0) return;
    event.preventDefault();
    event.stopPropagation();
    this.select(point.slot);
    this.drag = {
      kind: "point", pointerId: event.pointerId, slot: point.slot,
    };
    this.svg.setPointerCapture(event.pointerId);
  }

  beginAll(event) {
    if (event.button !== 0) return;
    event.preventDefault();
    event.stopPropagation();
    this.drag = {
      kind: "all", pointerId: event.pointerId,
      start: this.eventPosition(event), points: this.options.points(),
    };
    this.svg.setPointerCapture(event.pointerId);
  }

  dragPoint(position) {
    const points = this.options.points();
    const point = points.find(item => item.slot === this.drag.slot);
    if (!point) return;
    const frequency = point.fixed ? point.x :
      this.constrainFrequency(point.slot, this.frequency(position.x));
    this.options.setPoint(point.slot, frequency, this.logSeconds(position.y));
    this.options.select(point.slot);
    this.paint();
  }

  dragAll(position) {
    const delta = this.logSeconds(position.y) -
      this.logSeconds(this.drag.start.y);
    const points = this.drag.points.map(point => ({
      ...point,
      y: clamp(point.y + delta, this.options.minimumLogSeconds,
        this.options.maximumLogSeconds),
    }));
    this.options.replace(points, "body_decay_shift");
    this.paint();
  }

  select(slot) {
    this.selected = slot;
    this.options.select(slot);
    this.paint();
  }

  refresh() {
    const points = this.options.points();
    if (!points.some(point => point.slot === this.selected)) this.selected = 0;
    this.options.select(this.selected);
    this.paint();
  }

  paint() {
    this.svg.replaceChildren();
    this.paintGrid();
    const points = this.options.points();
    this.svg.append(element("polyline", {
      points: points.map(point =>
        `${this.xPosition(point.x)},${this.yPosition(point.y)}`).join(" "),
      class: "editor-curve decay-curve",
    }));
    for (const point of points) this.paintPoint(point);
    this.paintAllHandle(points);
    this.options.readout?.(`${points.length}/8 knots`);
  }

  paintGrid() {
    const nyquist = this.nyquist();
    const frequencies = [0, 100, 300, 1000, 3000, 10000, nyquist];
    for (const frequency of [...new Set(frequencies.filter(value =>
      value >= 0 && value <= nyquist))]) {
      const x = this.xPosition(frequency);
      this.svg.append(element("line", {
        x1: x, y1: View.top, x2: x, y2: View.height - View.bottom,
        class: "editor-grid",
      }));
      const label = element("text", {
        x, y: View.height - 11, class: "editor-tick", "text-anchor": "middle",
      });
      label.textContent = frequency === 0 ? "DC" : frequency === nyquist
        ? `Nyq ${Math.round(nyquist / 100) / 10}k`
        : frequency >= 1000 ? `${frequency / 1000}k` : frequency;
      this.svg.append(label);
    }
    for (const tick of this.options.yTicks) {
      const y = this.yPosition(tick.value);
      this.svg.append(element("line", {
        x1: View.left, y1: y, x2: View.width - View.right, y2: y,
        class: "editor-grid",
      }));
      const label = element("text", {
        x: View.left - 7, y: y + 3, class: "editor-tick",
        "text-anchor": "end",
      });
      label.textContent = tick.label;
      this.svg.append(label);
    }
  }

  paintPoint(point) {
    const node = element(point.fixed ? "rect" : "circle", point.fixed ? {
      x: this.xPosition(point.x) - 6, y: this.yPosition(point.y) - 6,
      width: 12, height: 12,
    } : {
      cx: this.xPosition(point.x), cy: this.yPosition(point.y), r: 6,
    });
    node.setAttribute("class",
      `editor-point decay-handle${point.slot === this.selected ? " selected" : ""}`);
    node.setAttribute("data-slot", point.slot);
    const tooltip = element("title");
    tooltip.textContent = point.fixed
      ? "Boundary frequency is fixed; drag vertically to set T60"
      : "Drag to move; double-click or press Delete to remove";
    node.append(tooltip);
    node.onpointerdown = event => this.beginPoint(event, point);
    node.ondblclick = event => {
      event.preventDefault(); event.stopPropagation();
      if (point.fixed) this.options.reset(point.slot);
      else this.options.remove(point.slot);
      this.refresh();
    };
    this.svg.append(node);
  }

  paintAllHandle(points) {
    const frequency = inverseErb(.5 * erb(this.nyquist()));
    const rate = erb(frequency);
    let right = 1;
    while (right < points.length && erb(points[right].x) < rate) ++right;
    right = Math.min(right, points.length - 1);
    const left = Math.max(0, right - 1);
    const span = erb(points[right].x) - erb(points[left].x);
    const amount = span > 1.e-6 ? (rate - erb(points[left].x)) / span : 0;
    const level = points[left].y + amount * (points[right].y - points[left].y);
    const x = this.xPosition(frequency); const y = this.yPosition(level);
    const handle = element("path", {
      d: `M ${x} ${y - 8} L ${x + 10} ${y} L ${x} ${y + 8} L ${x - 10} ${y} Z`,
      class: "decay-all-handle decay-handle",
      "aria-label": "Drag to move every T60 knot",
    });
    handle.onpointerdown = event => this.beginAll(event);
    this.svg.append(handle);
    const label = element("text", {
      x, y: y - 13, class: "decay-all-label", "text-anchor": "middle",
    });
    label.textContent = "ALL";
    this.svg.append(label);
  }
}
