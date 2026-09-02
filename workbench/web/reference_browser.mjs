import { readReferenceCorpora, readRemoteReference } from "./references.mjs";

const unique = values => [...new Set(values)];

export class ReferenceBrowser {
  constructor(elements, onReference, onStatus) {
    this.elements = elements;
    this.onReference = onReference;
    this.onStatus = onStatus;
    this.generation = 0;
  }

  async initialize() {
    this.corpora = await readReferenceCorpora();
    const { corpus, articulation, velocity, repeat } = this.elements;
    corpus.replaceChildren(...this.corpora.map(item => new Option(item.name, item.id)));
    if (!this.corpora.length) {
      corpus.add(new Option("No local corpus found", ""));
      return;
    }
    corpus.onchange = () => this.selectCorpus();
    articulation.onchange = () => {
      const velocityValue = this.velocity;
      const repeatValue = this.repeat;
      this.configureLayers(velocityValue, repeatValue);
      this.scheduleLoad();
    };
    velocity.oninput = () => { this.paint(); this.scheduleLoad(); };
    repeat.oninput = () => { this.paint(); this.scheduleLoad(); };
    corpus.ondblclick = event => {
      event.preventDefault();
      corpus.value = this.corpora[0].id;
      this.selectCorpus();
    };
    articulation.ondblclick = event => {
      event.preventDefault();
      const options = [...articulation.options].map(option => option.value);
      articulation.value = options.includes("edge") ? "edge" : options[0];
      this.configureLayers(96, 1);
      this.scheduleLoad();
    };
    velocity.ondblclick = event => {
      event.preventDefault();
      this.configureLayers(96, this.repeat);
      this.scheduleLoad();
    };
    repeat.ondblclick = event => {
      event.preventDefault();
      this.configureLayers(this.velocity, 1);
      this.scheduleLoad();
    };
    await this.selectCorpus();
  }

  async selectCorpus() {
    this.current = this.corpora.find(item => item.id === this.elements.corpus.value) ??
      this.corpora[0];
    const articulations = unique(this.current.cells.map(cell => cell.articulation));
    this.elements.articulation.replaceChildren(...articulations.map(
      value => new Option(value.replaceAll("-", " "), value),
    ));
    this.elements.articulation.value = articulations.includes("edge") ? "edge" : articulations[0];
    this.configureLayers(96, 1);
    await this.load();
  }

  configureLayers(preferredVelocity, preferredRepeat) {
    const cells = this.current.cells.filter(
      cell => cell.articulation === this.elements.articulation.value,
    );
    this.velocities = unique(cells.map(cell => cell.velocity)).sort((a, b) => a - b);
    this.repeats = unique(cells.map(cell => cell.repeat)).sort((a, b) => a - b);
    const velocityIndex = this.velocities.reduce((best, value, index) =>
      Math.abs(value - preferredVelocity) < Math.abs(this.velocities[best] - preferredVelocity)
        ? index : best, 0);
    Object.assign(this.elements.velocity, {
      min: 0, max: this.velocities.length - 1, step: 1, value: velocityIndex,
    });
    Object.assign(this.elements.repeat, {
      min: 0, max: this.repeats.length - 1, step: 1,
      value: Math.max(0, this.repeats.indexOf(preferredRepeat)),
    });
    this.paint();
  }

  paint() {
    this.elements.velocity.nextElementSibling.textContent =
      `v${String(this.velocity).padStart(3, "0")}`;
    this.elements.repeat.nextElementSibling.textContent = `take ${this.repeat}`;
  }

  scheduleLoad() {
    clearTimeout(this.timer);
    this.timer = setTimeout(() => this.load(), 100);
  }

  async load() {
    const generation = ++this.generation;
    const cell = this.current.cells.find(item =>
      item.articulation === this.elements.articulation.value &&
      item.velocity === this.velocity && item.repeat === this.repeat);
    if (!cell) return;
    this.onStatus(`Loading ${cell.label}…`);
    try {
      const reference = await readRemoteReference(this.current, cell);
      if (generation === this.generation) this.onReference(reference);
    } catch (error) {
      if (generation === this.generation) this.onStatus(String(error));
    }
  }

  get velocity() { return this.velocities[Number(this.elements.velocity.value)]; }
  get repeat() { return this.repeats[Number(this.elements.repeat.value)]; }
}
