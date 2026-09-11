import { svgPosition, decayDragDelta, shiftDecayPoints, decayPosition, decaySeconds } from "./decay_curve_geometry.mjs";
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
    this.width = View.width;
    this.resizeObserver = new ResizeObserver(() => {
      if (this.svg.isConnected) this.paint();
      else this.destroy();
    });
    this.resizeObserver.observe(this.svg);
    this.bind();
    this.paint();
  }

  destroy() { this.resizeObserver.disconnect(); }

  minimumFrequency() {
    return Math.max(0, Number(this.options.minimumFrequency));
  }

  maximumFrequency() {
    return Math.max(
      this.minimumFrequency() + 1, Number(this.options.maximumFrequency),
    );
  }

  xPosition(frequency) {
    const minimum = erb(this.minimumFrequency());
    const maximum = erb(this.maximumFrequency());
    const amount = (erb(clamp(
      frequency, this.minimumFrequency(), this.maximumFrequency(),
    )) - minimum) / (maximum - minimum);
    return View.left + amount * (this.width - View.left - View.right);
  }

  frequency(position) {
    const amount = clamp(
      (position - View.left) / (this.width - View.left - View.right), 0, 1,
    );
    const minimum = erb(this.minimumFrequency());
    const maximum = erb(this.maximumFrequency());
    return inverseErb(minimum + amount * (maximum - minimum));
  }

  yPosition(logSeconds) {
    const amount = decayPosition(2 ** logSeconds,
      2 ** this.options.minimumLogSeconds, 2 ** this.options.maximumLogSeconds);
    return View.height - View.bottom - amount *
      (View.height - View.top - View.bottom);
  }

  logSeconds(position) {
    const amount = clamp(
      (View.height - View.bottom - position) /
        (View.height - View.top - View.bottom), 0, 1,
    );
    return Math.log2(decaySeconds(amount,
      2 ** this.options.minimumLogSeconds, 2 ** this.options.maximumLogSeconds));
  }

  eventPosition(event) {
    return svgPosition(this.svg.getScreenCTM(), event.clientX, event.clientY);
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
      if (!position) return;
      if (this.drag.kind === "all") this.dragAll(position, event.shiftKey);
      else this.dragPoint(position, event.shiftKey);
      this.drag.start = position;
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
      if (!position) return;
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
    const start = this.eventPosition(event);
    if (!start) return;
    event.preventDefault();
    event.stopPropagation();
    this.select(point.slot);
    this.drag = {
      kind: "point", pointerId: event.pointerId, slot: point.slot, start,
    };
    this.svg.setPointerCapture(event.pointerId);
  }

  beginAll(event) {
    if (event.button !== 0) return;
    const start = this.eventPosition(event);
    if (!start) return;
    event.preventDefault();
    event.stopPropagation();
    this.drag = {
      kind: "all", pointerId: event.pointerId,
      start,
    };
    this.svg.setPointerCapture(event.pointerId);
  }

  dragPoint(position, fine = false) {
    const points = this.options.points();
    const point = points.find(item => item.slot === this.drag.slot);
    if (!point) return;
    const x = this.xPosition(point.x) + (position.x - this.drag.start.x) * (fine ? .1 : 1);
    const frequency = point.fixed ? point.x :
      this.constrainFrequency(point.slot, this.frequency(x));
    this.options.setPoint(point.slot, frequency, clamp(point.y + this.dragDelta(position, fine, point.y),
      this.options.minimumLogSeconds, this.options.maximumLogSeconds));
    this.options.select(point.slot);
    this.paint();
  }

  dragDelta(position, fine, anchor) {
    return decayDragDelta(this.drag.start.y, position.y,
      this.options.minimumLogSeconds, this.options.maximumLogSeconds,
      View.height - View.top - View.bottom, fine, anchor);
  }

  dragAll(position, fine = false) {
    const original = this.options.points();
    const points = shiftDecayPoints(original, this.dragDelta(position, fine, this.middleLevel(original)),
      this.options.minimumLogSeconds, this.options.maximumLogSeconds);
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
    this.width = Math.max(180, this.svg.clientWidth || View.width);
    this.svg.setAttribute("viewBox", `0 0 ${this.width} ${View.height}`);
    this.svg.replaceChildren();
    this.paintGrid();
    const points = this.options.points();
    this.svg.append(element("polyline", {
      points: this.curveCoordinates(points),
      class: "editor-curve decay-curve",
    }));
    for (const point of points) this.paintPoint(point);
    this.paintAllHandle(points);
    this.options.readout?.(`${points.length}/8 knots`);
  }

  // Sample the real ERB/log-T60 interpolation under the nonlinear display scale.
  curveCoordinates(points) {
    return points.slice(1).flatMap((right, index) => {
      const left = points[index], x = this.xPosition(left.x);
      const span = this.xPosition(right.x) - x;
      const steps = Math.max(1, Math.ceil(span / 3));
      return Array.from({ length: steps + 1 }, (_, step) => {
        const fraction = step / steps;
        return `${x + span * fraction},${this.yPosition(left.y + fraction * (right.y - left.y))}`;
      });
    }).join(" ");
  }

  paintGrid() {
    const minimum = this.minimumFrequency();
    const maximum = this.maximumFrequency();
    const frequencies = [minimum, 100, 300, 1000, 3000, 10000, maximum];
    for (const frequency of [...new Set(frequencies.filter(value =>
      value >= minimum && value <= maximum))]) {
      const x = this.xPosition(frequency);
      this.svg.append(element("line", {
        x1: x, y1: View.top, x2: x, y2: View.height - View.bottom,
        class: "editor-grid",
      }));
      const label = element("text", {
        x, y: View.height - 11, class: "editor-tick", "text-anchor": "middle",
      });
      label.textContent = frequency === minimum ? `${minimum} Hz` :
        frequency === maximum ? `${maximum / 1000}k`
        : frequency >= 1000 ? `${frequency / 1000}k` : frequency;
      this.svg.append(label);
    }
    for (const tick of this.options.yTicks) {
      const y = this.yPosition(tick.value);
      this.svg.append(element("line", {
        x1: View.left, y1: y, x2: this.width - View.right, y2: y,
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
      ? "Boundary frequency is fixed; drag vertically to set T60. Shift: 10× finer"
      : "Drag to move; Shift: 10× finer. Double-click or press Delete to remove";
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

  middleLevel(points) {
    const rate = .5 * (erb(this.minimumFrequency()) + erb(this.maximumFrequency()));
    let right = 1;
    while (right < points.length && erb(points[right].x) < rate) ++right;
    right = Math.min(right, points.length - 1);
    const left = Math.max(0, right - 1);
    const span = erb(points[right].x) - erb(points[left].x);
    const amount = span > 1.e-6 ? (rate - erb(points[left].x)) / span : 0;
    return points[left].y + amount * (points[right].y - points[left].y);
  }

  paintAllHandle(points) {
    const x = .5 * (View.left + this.width - View.right);
    const y = this.yPosition(this.middleLevel(points));
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
