const Svg = "http://www.w3.org/2000/svg";
const Height = 300;
const Plot = { left: 52, right: 18, top: 18, bottom: 32 };
const clamp = (value, minimum, maximum) =>
  Math.max(minimum, Math.min(maximum, value));

function erb(frequency) {
  return 21.4 * Math.log10(1 + .00437 * frequency);
}

function inverseErb(rate) {
  return (10 ** (rate / 21.4) - 1) / .00437;
}

function element(name, attributes = {}) {
  const result = document.createElementNS(Svg, name);
  for (const [key, value] of Object.entries(attributes)) {
    result.setAttribute(key, value);
  }
  return result;
}

export class ModalEditor {
  constructor(parent, options) {
    this.options = options;
    this.selected = null;
    this.drag = null;
    this.hover = null;
    this.tool = "edit";
    this.brushErb = 1;
    this.harmonicGuide = {
      visible: false, fundamentalHz: 110, snap: false,
    };
    this.width = 900;
    this.unitsPerPixel = 1;
    this.svg = element("svg", {
      class: "modal-editor", viewBox: `0 0 ${this.width} ${Height}`,
      role: "application", tabindex: 0,
      "aria-label": "Modal packet editor",
    });
    parent.append(this.svg);
    this.resizeObserver = new ResizeObserver(() => this.paint());
    this.resizeObserver.observe(parent);
    requestAnimationFrame(() => this.paint());
    this.bind();
    this.paint();
  }

  setTool(tool) {
    this.tool = ["edit", "shape", "paint"].includes(tool) ? tool : "edit";
    this.svg.classList.toggle("drawing", this.tool !== "edit");
    this.paint();
  }

  setBrushWidth(widthErb) {
    this.brushErb = clamp(Number(widthErb), .15, 5);
    this.paint();
  }

  setHarmonicGuide({ visible, fundamentalHz, snap }) {
    this.harmonicGuide = {
      visible: Boolean(visible),
      fundamentalHz: clamp(Number(fundamentalHz), 8, 8000),
      snap: Boolean(snap),
    };
    this.paint();
  }

  snapActiveToHarmonics() {
    if (!this.harmonicGuide.visible) return;
    const points = this.options.points();
    const fundamental = this.harmonicGuide.fundamentalHz;
    const first = Math.max(
      1, Math.ceil(this.options.minimumFrequency / fundamental),
    );
    const last = Math.floor(this.options.maximumFrequency / fundamental);
    const available = [];
    for (let harmonic = first; harmonic <= last; ++harmonic) {
      available.push(harmonic);
    }
    const active = points.map((point, index) => ({ point, index }))
      .filter(item => item.point.active)
      .sort((a, b) => a.point.frequency - b.point.frequency);
    for (const item of active) {
      if (!available.length) {
        item.point.frequency = this.nearestHarmonic(item.point.frequency);
        continue;
      }
      const target = item.point.frequency / fundamental;
      let choice = 0;
      for (let index = 1; index < available.length; ++index) {
        if (Math.abs(available[index] - target) <
            Math.abs(available[choice] - target)) choice = index;
      }
      item.point.frequency = available[choice] * fundamental;
      available.splice(choice, 1);
    }
    this.options.replace(points, "modal_harmonic_snap_all");
    this.options.select(this.selected);
    this.paint();
  }

  refresh() {
    if (this.selected !== null && !this.point(this.selected)?.active) {
      this.selected = null;
      this.options.select(null);
    }
    this.paint();
  }

  point(index) {
    return this.options.points()[index];
  }

  x(frequency) {
    frequency = clamp(
      frequency, this.options.minimumFrequency, this.options.maximumFrequency,
    );
    const first = erb(this.options.minimumFrequency);
    const last = erb(this.options.maximumFrequency);
    return Plot.left + (this.width - Plot.left - Plot.right) *
      (erb(frequency) - first) / (last - first);
  }

  frequency(position) {
    const amount = clamp(
      (position - Plot.left) / (this.width - Plot.left - Plot.right), 0, 1,
    );
    return inverseErb(
      erb(this.options.minimumFrequency) + amount *
      (erb(this.options.maximumFrequency) - erb(this.options.minimumFrequency)),
    );
  }

  snapFrequency(frequency) {
    frequency = clamp(
      frequency, this.options.minimumFrequency, this.options.maximumFrequency,
    );
    if (!this.harmonicGuide.visible || !this.harmonicGuide.snap) {
      return frequency;
    }
    return this.nearestHarmonic(frequency);
  }

  nearestHarmonic(frequency) {
    const fundamental = this.harmonicGuide.fundamentalHz;
    const first = Math.max(
      1, Math.ceil(this.options.minimumFrequency / fundamental),
    );
    const last = Math.floor(this.options.maximumFrequency / fundamental);
    if (last < first) return clamp(
      frequency, this.options.minimumFrequency, this.options.maximumFrequency,
    );
    return clamp(Math.round(frequency / fundamental), first, last) * fundamental;
  }

  y(level) {
    return Height - Plot.bottom - (Height - Plot.top - Plot.bottom) *
      (level - this.options.minimumLevel) /
      (this.options.maximumLevel - this.options.minimumLevel);
  }

  level(position) {
    const amount = clamp(
      (Height - Plot.bottom - position) /
      (Height - Plot.top - Plot.bottom), 0, 1,
    );
    return this.options.minimumLevel + amount *
      (this.options.maximumLevel - this.options.minimumLevel);
  }

  eventPosition(event) {
    const bounds = this.svg.getBoundingClientRect();
    return {
      x: this.width * (event.clientX - bounds.left) / bounds.width,
      y: Height * (event.clientY - bounds.top) / bounds.height,
    };
  }

  effectiveSpread(point) {
    const turbulence = clamp(
      this.options.globalTurbulence() * point.turbulence, 0, 1,
    );
    return turbulence * this.options.packetSpread();
  }

  bind() {
    this.svg.addEventListener("pointermove", event => {
      this.hover = { ...this.eventPosition(event), shift: event.shiftKey,
        control: event.ctrlKey || event.metaKey };
      if (this.drag?.pointerId === event.pointerId) this.updateDrag(event);
      else this.paint();
    });
    this.svg.addEventListener("pointerleave", () => {
      if (!this.drag) { this.hover = null; this.paint(); }
    });
    this.svg.addEventListener("pointerdown", event => {
      if (event.button !== 0 || event.target.closest(".modal-handle")) return;
      if (this.tool === "edit") {
        this.selected = null;
        this.options.select(null);
        this.paint();
        return;
      }
      event.preventDefault();
      this.drag = {
        kind: this.tool, pointerId: event.pointerId, lastPaint: null,
      };
      this.svg.setPointerCapture(event.pointerId);
      this.updateDrag(event);
    });
    this.svg.addEventListener("dblclick", event => {
      if (event.target.closest(".modal-handle")) return;
      event.preventDefault();
      const position = this.eventPosition(event);
      const index = this.options.insert(
        this.snapFrequency(this.frequency(position.x)), this.level(position.y),
      );
      if (index !== null) this.select(index);
    });
    const finish = event => {
      if (!this.drag || event.pointerId !== this.drag.pointerId) return;
      if (this.svg.hasPointerCapture(event.pointerId)) {
        this.svg.releasePointerCapture(event.pointerId);
      }
      this.drag = null;
      this.paint();
    };
    this.svg.addEventListener("pointerup", finish);
    this.svg.addEventListener("pointercancel", finish);
    this.svg.addEventListener("wheel", event => {
      if (this.selected === null) return;
      event.preventDefault();
      const points = this.options.points();
      const point = points[this.selected];
      point.turbulence = clamp(
        point.turbulence * Math.exp(-event.deltaY * .002), 0, 2,
      );
      this.options.replace(points, "modal_width");
      this.paint();
      this.options.select(this.selected);
    }, { passive: false });
    this.svg.addEventListener("keydown", event => {
      if ((event.key === "Delete" || event.key === "Backspace") &&
          this.selected !== null) {
        event.preventDefault();
        this.options.remove(this.selected);
        this.selected = null;
        this.options.select(null);
        this.paint();
      }
    });
  }

  beginHandle(event, index, kind) {
    if (event.button !== 0) return;
    event.preventDefault();
    event.stopPropagation();
    if (this.tool !== "edit") {
      this.drag = {
        kind: this.tool, pointerId: event.pointerId, lastPaint: null,
      };
      this.svg.setPointerCapture(event.pointerId);
      this.updateDrag(event);
      return;
    }
    this.select(index);
    const position = this.eventPosition(event);
    const point = this.point(index);
    this.drag = {
      kind: (event.ctrlKey || event.metaKey) ? "width" : kind,
      pointerId: event.pointerId, index, start: position,
      point: { ...point },
    };
    this.svg.setPointerCapture(event.pointerId);
  }

  updateDrag(event) {
    const position = this.eventPosition(event);
    if (this.drag.kind === "shape" || this.drag.kind === "paint") {
      const coalesced = event.getCoalescedEvents?.();
      for (const sample of coalesced?.length ? coalesced : [event]) {
        this.paintAlong(this.eventPosition(sample), sample);
      }
      return;
    }
    const points = this.options.points();
    const point = points[this.drag.index];
    if (!point) return;
    if (this.drag.kind === "width") {
      const distance = Math.abs(erb(this.frequency(position.x)) -
        erb(this.drag.point.frequency));
      const global = Math.max(
        .05, this.options.globalTurbulence() * this.options.packetSpread(),
      );
      point.turbulence = clamp(distance / global, 0, 2);
    } else {
      const precision = event.shiftKey ? .2 : 1;
      point.frequency = this.snapFrequency(clamp(
        this.drag.point.frequency + precision *
        (this.frequency(position.x) - this.drag.point.frequency),
        this.options.minimumFrequency, this.options.maximumFrequency,
      ));
      point.level = clamp(
        this.drag.point.level + precision *
        (this.level(position.y) - this.drag.point.level),
        this.options.minimumLevel, this.options.maximumLevel,
      );
    }
    this.options.replace(points, this.drag.kind === "width"
      ? "modal_width" : "resolved_modes");
    this.options.select(this.drag.index);
    this.paint();
  }

  paintAlong(position, event) {
    const previous = this.drag.lastPaint;
    const samples = [];
    if (!previous) {
      samples.push(position);
    } else {
      const firstErb = erb(this.frequency(previous.x));
      const lastErb = erb(this.frequency(position.x));
      const spacing = clamp(.65 * this.brushErb, .3, 1.25);
      const steps = Math.max(
        1, Math.ceil(Math.abs(lastErb - firstErb) / spacing),
      );
      for (let step = 1; step <= steps; ++step) {
        const amount = step / steps;
        samples.push({
          x: previous.x + amount * (position.x - previous.x),
          y: previous.y + amount * (position.y - previous.y),
        });
      }
    }
    const points = this.options.points();
    for (const sample of samples) this.paintAt(sample, event, points);
    this.drag.lastPaint = position;
    this.options.replace(points, "resolved_modes");
    if (this.selected !== null) this.options.select(this.selected);
    this.paint();
  }

  paintAt(position, event, points) {
    const centre = erb(this.frequency(position.x));
    const modifier = event.ctrlKey || event.metaKey ? 2.5 :
      event.shiftKey ? .35 : 1;
    const sigma = this.brushErb * modifier;
    const target = this.level(position.y);
    if (this.drag.kind === "paint") {
      const candidateFrequency = this.snapFrequency(this.frequency(position.x));
      const candidateErb = erb(candidateFrequency);
      const placementRadius = clamp(.7 * sigma, .3, 1.4);
      const occupied = points.some(point => point.active &&
        Math.abs(erb(point.frequency) - candidateErb) < placementRadius);
      const index = occupied ? -1 : points.findIndex(point => !point.active);
      if (index >= 0) {
        points[index] = {
          frequency: candidateFrequency,
          level: Math.max(target, this.options.minimumLevel + .1),
          turbulence: 1,
          active: true,
        };
        this.selected = index;
      }
    }
    for (const point of points) {
      if (!point.active) continue;
      const distance = (erb(point.frequency) - centre) / sigma;
      const weight = Math.exp(-.5 * distance * distance);
      point.level = clamp(
        point.level + .65 * weight * (target - point.level),
        this.options.minimumLevel, this.options.maximumLevel,
      );
    }
  }

  select(index) {
    this.selected = index;
    this.options.select(index);
    this.paint();
  }

  syncGeometry() {
    const bounds = this.svg.getBoundingClientRect();
    if (!bounds.width || !bounds.height) return;
    this.unitsPerPixel = Height / bounds.height;
    const width = Math.max(600, Height * bounds.width / bounds.height);
    if (Math.abs(width - this.width) < 1) return;
    this.width = width;
    this.svg.setAttribute("viewBox", `0 0 ${this.width} ${Height}`);
  }

  paint() {
    this.syncGeometry();
    this.svg.replaceChildren();
    this.paintGrid();
    this.paintHarmonicGrid();
    const points = this.options.points();
    points.forEach((point, index) => {
      if (point.active) this.paintPacket(point, index);
    });
    if (this.hover && this.tool !== "edit") this.paintBrushPreview();
    const count = points.filter(point => point.active).length;
    this.options.readout(
      this.selected === null ? `${count}/${points.length} active modes` :
        this.pointText(points[this.selected]),
    );
  }

  paintGrid() {
    for (const frequency of [50, 100, 300, 1000, 3000, 6000, 10000, 15000]) {
      const x = this.x(frequency);
      this.svg.append(element("line", {
        x1: x, y1: Plot.top, x2: x, y2: Height - Plot.bottom,
        class: "editor-grid",
      }));
      const label = element("text", {
        x, y: Height - 10, class: "editor-tick", "text-anchor": "middle",
      });
      label.style.fontSize = `${9 * this.unitsPerPixel}px`;
      label.textContent = frequency >= 1000 ? `${frequency / 1000}k` : frequency;
      this.svg.append(label);
    }
    for (const level of [-72, -48, -24, 0, 6]) {
      const y = this.y(level);
      this.svg.append(element("line", {
        x1: Plot.left, y1: y, x2: this.width - Plot.right, y2: y,
        class: "editor-grid",
      }));
      const label = element("text", {
        x: Plot.left - 7, y: y + 3, class: "editor-tick",
        "text-anchor": "end",
      });
      label.style.fontSize = `${9 * this.unitsPerPixel}px`;
      label.textContent = level === -72 ? "off" : `${level} dB`;
      this.svg.append(label);
    }
  }

  paintHarmonicGrid() {
    if (!this.harmonicGuide.visible) return;
    const fundamental = this.harmonicGuide.fundamentalHz;
    const first = Math.max(
      1, Math.ceil(this.options.minimumFrequency / fundamental),
    );
    const last = Math.floor(this.options.maximumFrequency / fundamental);
    const lineGap = 7 * this.unitsPerPixel;
    const labelGap = 42 * this.unitsPerPixel;
    let previousLine = -Infinity;
    let previousLabel = -Infinity;
    for (let harmonic = first; harmonic <= last; ++harmonic) {
      const frequency = harmonic * fundamental;
      const x = this.x(frequency);
      if (x - previousLine < lineGap && harmonic !== 1) continue;
      const primary = harmonic === 1;
      this.svg.append(element("line", {
        x1: x, y1: Plot.top, x2: x, y2: Height - Plot.bottom,
        class: `harmonic-grid${primary ? " fundamental" : ""}`,
        "data-harmonic": harmonic,
      }));
      previousLine = x;
      if (x - previousLabel < labelGap && !primary) continue;
      const label = element("text", {
        x: x + 3 * this.unitsPerPixel,
        y: Plot.top + 11 * this.unitsPerPixel,
        class: "harmonic-label", "data-harmonic": harmonic,
      });
      label.style.fontSize = `${9 * this.unitsPerPixel}px`;
      label.textContent = `${harmonic}×`;
      this.svg.append(label);
      previousLabel = x;
    }
  }

  paintPacket(point, index) {
    const centreErb = erb(point.frequency);
    const spread = this.effectiveSpread(point);
    if (spread > .025) {
      const sigma = Math.max(.08, .45 * spread);
      const samples = [];
      for (let step = -24; step <= 24; ++step) {
        const offset = 3 * sigma * step / 24;
        const envelope = Math.exp(-.5 * (offset / sigma) ** 2);
        const level = this.options.minimumLevel + envelope *
          (point.level - this.options.minimumLevel);
        samples.push(`${this.x(inverseErb(centreErb + offset))},${this.y(level)}`);
      }
      const path = `${this.x(inverseErb(centreErb - 3 * sigma))},${this.y(this.options.minimumLevel)} ` +
        samples.join(" ") + ` ${this.x(inverseErb(centreErb + 3 * sigma))},${this.y(this.options.minimumLevel)}`;
      this.svg.append(element("polygon", {
        points: path,
        class: `modal-packet${index === this.selected ? " selected" : ""}`,
      }));
    }
    const bar = element("line", {
      x1: this.x(point.frequency), y1: this.y(this.options.minimumLevel),
      x2: this.x(point.frequency), y2: this.y(point.level),
      class: `modal-bar modal-handle${index === this.selected ? " selected" : ""}`,
    });
    const node = element("circle", {
      cx: this.x(point.frequency), cy: this.y(point.level),
      r: 6 * this.unitsPerPixel,
      class: `modal-node modal-handle${index === this.selected ? " selected" : ""}`,
    });
    for (const handle of [bar, node]) {
      handle.onpointerdown = event => this.beginHandle(event, index, "centre");
      handle.ondblclick = event => {
        event.preventDefault(); event.stopPropagation();
        this.options.remove(index); this.selected = null;
        this.options.select(null); this.paint();
      };
    }
    this.svg.append(bar, node);
    if (index === this.selected && spread > .025) {
      for (const direction of [-1, 1]) {
        const wing = element("circle", {
          cx: this.x(inverseErb(centreErb + direction * spread)),
          cy: this.y(this.options.minimumLevel +
            .45 * (point.level - this.options.minimumLevel)),
          r: 5 * this.unitsPerPixel, class: "modal-wing modal-handle",
        });
        wing.onpointerdown = event => this.beginHandle(event, index, "width");
        this.svg.append(wing);
      }
    }
  }

  paintBrushPreview() {
    const centre = erb(this.frequency(this.hover.x));
    const modifier = this.hover.control ? 2.5 : this.hover.shift ? .35 : 1;
    const sigma = this.brushErb * modifier;
    const target = this.level(this.hover.y);
    const samples = [];
    for (let step = -24; step <= 24; ++step) {
      const offset = 3 * sigma * step / 24;
      const envelope = Math.exp(-.5 * (offset / sigma) ** 2);
      const level = this.options.minimumLevel + envelope *
        (target - this.options.minimumLevel);
      samples.push(`${this.x(inverseErb(centre + offset))},${this.y(level)}`);
    }
    this.svg.append(element("polyline", {
      points: samples.join(" "),
      class: `modal-brush-preview ${this.tool}`,
    }));
  }

  pointText(point) {
    return `${point.frequency < 1000 ? point.frequency.toFixed(1) :
      (point.frequency / 1000).toFixed(3) + " k"}Hz · ` +
      `${point.level.toFixed(1)} dB · turbulence ${point.turbulence.toFixed(2)}x`;
  }
}
