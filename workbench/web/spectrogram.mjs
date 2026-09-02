const magma = [
  [0, 0, 4], [28, 16, 68], [79, 18, 123], [129, 37, 129],
  [181, 54, 122], [229, 80, 100], [251, 135, 97], [254, 194, 135],
  [252, 253, 191],
];

function colour(stops, amount) {
  const position = Math.max(0, Math.min(1, amount)) * (stops.length - 1);
  const first = Math.floor(position);
  const second = Math.min(stops.length - 1, first + 1);
  const mix = position - first;
  return stops[first].map((value, index) =>
    Math.round(value + mix * (stops[second][index] - value)));
}

function differenceColour(value, span) {
  const amount = Math.max(-1, Math.min(1, value / span));
  return amount < 0
    ? colour([[16, 48, 112], [238, 242, 245]], amount + 1)
    : colour([[238, 242, 245], [154, 24, 54]], amount);
}

function timeLabel(seconds, span) {
  if (Math.abs(seconds) < 0.0005) seconds = 0;
  if (span < 0.1) return `${(1000 * seconds).toFixed(0)} ms`;
  if (span < 2) return `${seconds.toFixed(2)} s`;
  return `${seconds.toFixed(1)} s`;
}

function buildPalette(colourAt, size = 256) {
  const result = new Uint8ClampedArray(size * 3);
  for (let index = 0; index < size; ++index) {
    const rgb = colourAt(index / (size - 1));
    result.set(rgb, index * 3);
  }
  return result;
}

const magmaPalette = buildPalette(amount => colour(magma, amount));
const differencePalette = buildPalette(amount =>
  differenceColour(2 * amount - 1, 1));

export function wheelPanSeconds(
  deltaX, deltaY, deltaMode, width, height, timeSpan,
) {
  const scale = deltaMode === 1 ? 16 : deltaMode === 2 ? height : 1;
  const primary = Math.abs(deltaY) >= Math.abs(deltaX) ? deltaY : deltaX;
  return primary * scale / Math.max(1, width) * timeSpan * 0.25;
}

export class SpectrogramView {
  constructor(canvas) {
    this.canvas = canvas;
    this.context = canvas.getContext("2d", { alpha: false });
    this.reference = null;
    this.synthesis = null;
    this.settings = {
      mode: "mirror", floorDb: -160, dynamicRangeDb: 90,
      differenceDb: 24, split: 0.5,
    };
    this.viewport = { start: 0, end: 1, lowHz: 40, highHz: 24000 };
    this.timeOffsets = { reference: 0, synthesis: 0 };
    canvas.title = "Mirror: wheel pans both sides; drag one side to align; Ctrl+wheel zooms; double-click resets";
    this.#bindInteraction();
    new ResizeObserver(() => this.draw()).observe(canvas);
  }

  setData(kind, data) {
    this[kind] = data;
    const duration = data.frames * data.hop / data.sampleRate;
    this.viewport.end = Math.max(this.viewport.end, duration);
    this.viewport.highHz = Math.min(this.viewport.highHz, data.sampleRate / 2);
    this.draw();
  }

  setSettings(settings) { Object.assign(this.settings, settings); this.draw(); }
  reset() {
    const data = this.reference ?? this.synthesis;
    if (!data) return;
    this.viewport = {
      start: 0, end: data.frames * data.hop / data.sampleRate,
      lowHz: 40, highHz: data.sampleRate / 2,
    };
    this.timeOffsets.reference = 0;
    this.timeOffsets.synthesis = 0;
    this.draw();
    this.onViewport?.(this.viewport);
  }

  #frameMap(data, width, offset = 0) {
    const result = new Int32Array(width); result.fill(-1);
    if (!data) return result;
    const duration = (data.frames - 1) * data.hop / data.sampleRate;
    for (let x = 0; x < width; ++x) {
      const time = this.viewport.start + offset + x / width *
        (this.viewport.end - this.viewport.start);
      if (time >= 0 && time <= duration) result[x] = Math.max(
        0, Math.min(data.frames - 1, Math.round(time * data.sampleRate / data.hop)),
      );
    }
    return result;
  }

  #binMap(data, height) {
    const result = new Int32Array(height); result.fill(-1);
    if (!data) return result;
    const low = Math.log(Math.max(1, this.viewport.lowHz));
    const high = Math.log(Math.max(2, this.viewport.highHz));
    for (let y = 0; y < height; ++y) {
      const frequency = Math.exp(low + (1 - y / height) * (high - low));
      result[y] = Math.max(0, Math.min(data.bins - 1,
        Math.round(frequency * data.size / data.sampleRate)));
    }
    return result;
  }

  #source(mode, x, y, width, height) {
    if (mode === "reference") return [this.reference, x, y];
    if (mode === "synthesis") return [this.synthesis, x, y];
    if (mode === "horizontal") {
      const divider = height * this.settings.split;
      const top = y < divider;
      const local = top ? y / divider : (y - divider) / (height - divider);
      return [top ? this.reference : this.synthesis, x, local * height];
    }
    if (mode === "vertical") {
      const divider = width * this.settings.split;
      const left = x < divider;
      const local = left ? x / divider : (x - divider) / (width - divider);
      return [left ? this.reference : this.synthesis, local * width, y];
    }
    const divider = width * this.settings.split;
    const left = x < divider;
    const mirroredX = left ? (divider - x - 1) / divider * width :
      (x - divider) / (width - divider) * width;
    return [left ? this.reference : this.synthesis, mirroredX, y];
  }

  #divider(width, height) {
    if (this.settings.mode === "vertical" || this.settings.mode === "mirror") {
      return { axis: "x", position: width * this.settings.split };
    }
    if (this.settings.mode === "horizontal") {
      return { axis: "y", position: height * this.settings.split };
    }
    return null;
  }

  #value(data, frames, bins, x, y, width, height) {
    if (!data) return this.settings.floorDb;
    const frame = frames[Math.max(0, Math.min(width - 1, Math.floor(x)))];
    const bin = bins[Math.max(0, Math.min(height - 1, Math.floor(y)))];
    return frame < 0 || bin < 0
      ? this.settings.floorDb : data.values[frame * data.bins + bin];
  }

  #drawLegend(width, ratio, divider) {
    const entries = this.settings.mode === "reference"
      ? [["Reference", "#e8b45a"]]
      : this.settings.mode === "synthesis"
        ? [["TriggerFish", "#68a7d8"]]
        : this.settings.mode === "difference"
          ? [["TriggerFish − reference", "#e7e9ed"]]
          : [["Reference", "#e8b45a"], ["TriggerFish", "#68a7d8"]];
    const context = this.context;
    context.font = `${11 * ratio}px sans-serif`;
    context.textBaseline = "top";
    const gap = 18 * ratio;
    const widths = entries.map(([label]) => context.measureText(label).width + gap);
    const total = widths.reduce((sum, value) => sum + value, 0) +
      (entries.length - 1) * 14 * ratio;
    let x = Math.max(8 * ratio, (width - total) / 2);
    context.fillStyle = "rgba(8, 11, 15, 0.78)";
    context.fillRect(x - 6 * ratio, 5 * ratio,
                     total + 12 * ratio, 20 * ratio);
    entries.forEach(([label, colourValue], index) => {
      context.fillStyle = colourValue;
      context.fillRect(x, 10 * ratio, 8 * ratio, 8 * ratio);
      context.fillStyle = "rgba(240, 243, 247, 0.94)";
      context.fillText(label, x + 12 * ratio, 7 * ratio);
      x += widths[index] + 14 * ratio;
    });

    if (this.settings.mode !== "mirror" || !divider) return;
    const signed = value => `${value >= 0 ? "+" : ""}${value.toFixed(3)} s`;
    context.fillStyle = "rgba(235, 240, 246, 0.82)";
    context.textAlign = "right";
    context.fillText(`offset ${signed(this.timeOffsets.reference)}`,
                     divider.position - 7 * ratio, 27 * ratio);
    context.textAlign = "left";
    context.fillText(`offset ${signed(this.timeOffsets.synthesis)}`,
                     divider.position + 7 * ratio, 27 * ratio);
    context.textAlign = "start";
  }

  #drawTimeAxis(width, height, ratio, divider) {
    const context = this.context;
    const axisHeight = 25 * ratio;
    const baseline = height - axisHeight;
    const span = this.viewport.end - this.viewport.start;
    context.fillStyle = "rgba(8, 11, 15, 0.82)";
    context.fillRect(0, baseline, width, axisHeight);
    context.strokeStyle = "rgba(235, 240, 246, 0.58)";
    context.fillStyle = "rgba(235, 240, 246, 0.86)";
    context.lineWidth = ratio;
    context.font = `${10 * ratio}px sans-serif`;
    context.textBaseline = "top";

    const tick = (x, time, alignment, labelInset = 0) => {
      context.beginPath();
      context.moveTo(Math.round(x), baseline);
      context.lineTo(Math.round(x), baseline + 5 * ratio);
      context.stroke();
      context.textAlign = alignment;
      context.fillText(
        timeLabel(time, span), x + labelInset, baseline + 7 * ratio,
      );
    };
    const divisions = 4;
    if (this.settings.mode === "mirror" && divider) {
      for (let index = 0; index <= divisions; ++index) {
        const amount = index / divisions;
        const referenceX = divider.position * (1 - amount);
        const synthesisX = divider.position +
          (width - divider.position) * amount;
        const referenceTime = this.viewport.start +
          this.timeOffsets.reference + amount * span;
        const synthesisTime = this.viewport.start +
          this.timeOffsets.synthesis + amount * span;
        tick(referenceX, referenceTime,
             index === 0 ? "right" : index === divisions ? "left" : "center",
             index === 0 ? -4 * ratio : 0);
        tick(synthesisX, synthesisTime,
             index === 0 ? "left" : index === divisions ? "right" : "center",
             index === 0 ? 4 * ratio : 0);
      }
    } else {
      for (let index = 0; index <= divisions; ++index) {
        const amount = index / divisions;
        tick(amount * width, this.viewport.start + amount * span,
             index === 0 ? "left" : index === divisions ? "right" : "center");
      }
    }
    context.textAlign = "start";
  }

  draw() {
    const ratio = devicePixelRatio || 1;
    const width = Math.max(1, Math.round(this.canvas.clientWidth * ratio));
    const height = Math.max(1, Math.round(this.canvas.clientHeight * ratio));
    if (this.canvas.width !== width || this.canvas.height !== height) {
      this.canvas.width = width; this.canvas.height = height;
    }
    const image = this.context.createImageData(width, height);
    const mode = this.settings.mode;
    const availablePeaks = [this.reference?.peakDb, this.synthesis?.peakDb]
      .filter(Number.isFinite);
    const ceiling = availablePeaks.length ? Math.max(...availablePeaks) : 0;
    const floor = ceiling - this.settings.dynamicRangeDb;
    const referenceFrames = this.#frameMap(
      this.reference, width, this.timeOffsets.reference,
    );
    const synthesisFrames = this.#frameMap(
      this.synthesis, width, this.timeOffsets.synthesis,
    );
    const referenceBins = this.#binMap(this.reference, height);
    const synthesisBins = this.#binMap(this.synthesis, height);
    for (let y = 0; y < height; ++y) {
      for (let x = 0; x < width; ++x) {
        let palette = magmaPalette;
        let colourIndex;
        if (mode === "difference") {
          const reference = this.#value(
            this.reference, referenceFrames, referenceBins, x, y, width, height,
          );
          const synthesis = this.#value(
            this.synthesis, synthesisFrames, synthesisBins, x, y, width, height,
          );
          palette = differencePalette;
          colourIndex = Math.round(255 * (0.5 + 0.5 * Math.max(-1, Math.min(
            1, (synthesis - reference) / this.settings.differenceDb,
          ))));
        } else {
          const source = this.#source(mode, x, y, width, height);
          const isReference = source[0] === this.reference;
          const value = this.#value(source[0],
            isReference ? referenceFrames : synthesisFrames,
            isReference ? referenceBins : synthesisBins,
            source[1], source[2], width, height);
          colourIndex = Math.round(255 * Math.max(0, Math.min(
            1, (value - floor) / (ceiling - floor),
          )));
        }
        const offset = 4 * (y * width + x);
        const paletteOffset = 3 * colourIndex;
        image.data[offset] = palette[paletteOffset];
        image.data[offset + 1] = palette[paletteOffset + 1];
        image.data[offset + 2] = palette[paletteOffset + 2];
        image.data[offset + 3] = 255;
      }
    }
    this.context.putImageData(image, 0, 0);
    const divider = this.#divider(width, height);
    if (divider) {
      this.context.fillStyle = "rgba(255, 255, 255, 0.72)";
      if (divider.axis === "x") {
        this.context.fillRect(Math.round(divider.position) - 1, 0, 2, height);
      } else {
        this.context.fillRect(0, Math.round(divider.position) - 1, width, 2);
      }
    }
    this.#drawLegend(width, ratio, divider);
    this.#drawTimeAxis(width, height, ratio, divider);
  }

  #bindInteraction() {
    let drag = null;
    this.canvas.addEventListener("wheel", event => {
      event.preventDefault();
      if (this.settings.mode === "mirror" && !event.ctrlKey && !event.metaKey) {
        const shift = wheelPanSeconds(
          event.deltaX, event.deltaY, event.deltaMode,
          this.canvas.clientWidth, this.canvas.clientHeight,
          this.viewport.end - this.viewport.start,
        );
        this.timeOffsets.reference += shift;
        this.timeOffsets.synthesis += shift;
        this.draw();
        return;
      }
      const amount = Math.exp(event.deltaY * 0.001);
      const position = event.offsetX / this.canvas.clientWidth;
      const span = (this.viewport.end - this.viewport.start) * amount;
      const centre = this.viewport.start + position *
        (this.viewport.end - this.viewport.start);
      this.viewport.start = centre - position * span;
      this.viewport.end = this.viewport.start + Math.max(0.002, span);
      this.draw();
      this.onViewport?.(this.viewport);
    }, { passive: false });
    this.canvas.addEventListener("pointerdown", event => {
      const divider = this.#divider(
        this.canvas.clientWidth, this.canvas.clientHeight,
      );
      const coordinate = divider?.axis === "x" ? event.offsetX : event.offsetY;
      if (divider && Math.abs(coordinate - divider.position) <= 8) {
        drag = { kind: "split", axis: divider.axis };
        this.canvas.setPointerCapture(event.pointerId);
        return;
      }
      if (this.settings.mode === "mirror") {
        const source = event.offsetX < divider.position
          ? "reference" : "synthesis";
        drag = {
          kind: "source-pan", source, x: event.clientX,
          offset: this.timeOffsets[source],
        };
        this.canvas.setPointerCapture(event.pointerId);
        return;
      }
      drag = { x: event.clientX, start: this.viewport.start, end: this.viewport.end };
      this.canvas.setPointerCapture(event.pointerId);
    });
    this.canvas.addEventListener("pointermove", event => {
      if (!drag) {
        const divider = this.#divider(
          this.canvas.clientWidth, this.canvas.clientHeight,
        );
        const coordinate = divider?.axis === "x" ? event.offsetX : event.offsetY;
        const near = divider && Math.abs(coordinate - divider.position) <= 8;
        this.canvas.style.cursor = near
          ? (divider.axis === "x" ? "col-resize" : "row-resize") : "grab";
        return;
      }
      if (drag.kind === "split") {
        const bounds = this.canvas.getBoundingClientRect();
        const amount = drag.axis === "x"
          ? (event.clientX - bounds.left) / bounds.width
          : (event.clientY - bounds.top) / bounds.height;
        this.settings.split = Math.max(0.1, Math.min(0.9, amount));
        this.draw();
        return;
      }
      if (drag.kind === "source-pan") {
        const divider = this.canvas.clientWidth * this.settings.split;
        const sideWidth = drag.source === "reference"
          ? divider : this.canvas.clientWidth - divider;
        const direction = drag.source === "reference" ? 1 : -1;
        const span = this.viewport.end - this.viewport.start;
        const shift = direction * (event.clientX - drag.x) /
          Math.max(1, sideWidth) * span;
        this.timeOffsets[drag.source] = drag.offset + shift;
        this.draw();
        return;
      }
      const span = drag.end - drag.start;
      const shift = -(event.clientX - drag.x) / this.canvas.clientWidth * span;
      this.viewport.start = Math.max(0, drag.start + shift);
      this.viewport.end = this.viewport.start + span;
      this.draw();
      this.onViewport?.(this.viewport);
    });
    const endDrag = () => { drag = null; };
    this.canvas.addEventListener("pointerup", endDrag);
    this.canvas.addEventListener("pointercancel", endDrag);
    this.canvas.addEventListener("dblclick", () => this.reset());
  }
}
