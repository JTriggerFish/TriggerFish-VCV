const Svg = "http://www.w3.org/2000/svg";
const clamp = (value, minimum, maximum) =>
  Math.max(minimum, Math.min(maximum, value));

function svgElement(name, attributes = {}) {
  const element = document.createElementNS(Svg, name);
  for (const [key, value] of Object.entries(attributes)) {
    element.setAttribute(key, value);
  }
  return element;
}

function erb(frequency) {
  return 21.4 * Math.log10(1 + .00437 * frequency);
}

function inverseErb(rate) {
  return (10 ** (rate / 21.4) - 1) / .00437;
}

export class PointEditor {
  constructor(parent, options) {
    this.options = options;
    this.selected = 0;
    this.svg = svgElement("svg", {
      class: "point-editor", viewBox: "0 0 360 180", role: "img",
      "aria-label": options.label,
    });
    parent.append(this.svg);
    this.bindPainting();
    this.paint();
  }

  bindPainting() {
    if (!this.options.paint) return;
    this.svg.oncontextmenu = event => event.preventDefault();
    this.svg.addEventListener("pointerdown", event => {
      if (event.target !== this.svg &&
          !event.target.classList.contains("editor-grid")) return;
      const erase = event.shiftKey || event.button === 2;
      const update = current => {
        const bounds = this.svg.getBoundingClientRect();
        const x = 360 * (current.clientX - bounds.left) / bounds.width;
        const y = 180 * (current.clientY - bounds.top) / bounds.height;
        const points = this.options.points();
        let nearest = 0;
        for (let index = 1; index < points.length; ++index) {
          if (Math.abs(this.xPosition(points[index].x) - x) <
              Math.abs(this.xPosition(points[nearest].x) - x)) nearest = index;
        }
        this.selected = nearest;
        this.options.select?.(nearest);
        this.options.setPoint(
          nearest, points[nearest].x,
          erase ? this.options.yMinimum : this.yValue(y),
        );
        this.paint();
      };
      this.svg.setPointerCapture(event.pointerId);
      update(event);
      this.svg.onpointermove = update;
    });
    this.svg.addEventListener("pointerup", () => {
      this.svg.onpointermove = null;
    });
  }

  xPosition(value) {
    const { xMinimum, xMaximum, xScale = "log" } = this.options;
    const first = xScale === "erb" ? erb(xMinimum) : Math.log(xMinimum);
    const last = xScale === "erb" ? erb(xMaximum) : Math.log(xMaximum);
    const current = xScale === "erb" ? erb(value) : Math.log(value);
    return 34 + 308 * (current - first) / (last - first);
  }

  xValue(position) {
    const amount = clamp((position - 34) / 308, 0, 1);
    const { xMinimum, xMaximum, xScale = "log" } = this.options;
    if (xScale === "erb") {
      return inverseErb(erb(xMinimum) + amount *
        (erb(xMaximum) - erb(xMinimum)));
    }
    return Math.exp(Math.log(xMinimum) + amount *
      Math.log(xMaximum / xMinimum));
  }

  yPosition(value) {
    const { yMinimum, yMaximum } = this.options;
    return 148 - 124 * (value - yMinimum) / (yMaximum - yMinimum);
  }

  yValue(position) {
    const amount = clamp((148 - position) / 124, 0, 1);
    return this.options.yMinimum + amount *
      (this.options.yMaximum - this.options.yMinimum);
  }

  paint() {
    this.svg.replaceChildren();
    const points = this.options.points();
    for (const frequency of [100, 300, 1000, 3000, 10000]) {
      if (frequency < this.options.xMinimum || frequency > this.options.xMaximum)
        continue;
      const x = this.xPosition(frequency);
      this.svg.append(svgElement("line", {
        x1: x, y1: 20, x2: x, y2: 150, class: "editor-grid",
      }));
      const label = svgElement("text", {
        x, y: 166, class: "editor-tick", "text-anchor": "middle",
      });
      label.textContent = frequency >= 1000 ? `${frequency / 1000}k` : frequency;
      this.svg.append(label);
    }
    for (const tick of this.options.yTicks) {
      const y = this.yPosition(tick.value);
      this.svg.append(svgElement("line", {
        x1: 32, y1: y, x2: 344, y2: y, class: "editor-grid",
      }));
      const label = svgElement("text", {
        x: 28, y: y + 3, class: "editor-tick", "text-anchor": "end",
      });
      label.textContent = tick.label;
      this.svg.append(label);
    }
    if (this.options.connected) {
      this.svg.append(svgElement("polyline", {
        points: points.map(point =>
          `${this.xPosition(point.x)},${this.yPosition(point.y)}`).join(" "),
        class: "editor-curve",
      }));
    }
    points.forEach((point, index) => this.addPoint(point, index));
  }

  addPoint(point, index) {
    const bindDrag = element => {
      element.onpointerdown = event => {
        event.stopPropagation();
        this.selected = index;
        this.options.select?.(index);
        this.svg.setPointerCapture(event.pointerId);
        update(event);
        this.svg.onpointermove = update;
      };
      element.ondblclick = event => {
        event.preventDefault(); event.stopPropagation();
        this.options.resetPoint(index);
        this.paint();
      };
    };
    const update = event => {
      const bounds = this.svg.getBoundingClientRect();
      const x = 360 * (event.clientX - bounds.left) / bounds.width;
      const y = 180 * (event.clientY - bounds.top) / bounds.height;
      const nextX = this.options.movableX ? this.xValue(x) : point.x;
      this.options.setPoint(index, nextX, this.yValue(y));
      this.paint();
    };
    if (this.options.bars) {
      const bar = svgElement("line", {
        x1: this.xPosition(point.x), y1: this.yPosition(this.options.yMinimum),
        x2: this.xPosition(point.x), y2: this.yPosition(point.y),
        class: `editor-bar${index === this.selected ? " selected" : ""}`,
      });
      bindDrag(bar);
      this.svg.append(bar);
    }
    const circle = svgElement("circle", {
      cx: this.xPosition(point.x), cy: this.yPosition(point.y), r: 5,
      class: `editor-point${index === this.selected ? " selected" : ""}`,
      tabindex: 0,
    });
    this.svg.onpointerup = () => { this.svg.onpointermove = null; };
    bindDrag(circle);
    this.svg.append(circle);
  }
}
