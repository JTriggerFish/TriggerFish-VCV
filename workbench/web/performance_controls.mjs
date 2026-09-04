const byId = id => document.getElementById(id);
const clamp01 = value => Math.max(0, Math.min(1, value));

const implementFamilies = [
  {
    value: 0, label: "Bristle stiffness",
    endpoints: ["soft / fine", "stiff / coarse"],
  },
  {
    value: .5, label: "Mallet firmness",
    endpoints: ["soft", "hard"],
  },
  {
    value: 1, label: "Tip hardness",
    endpoints: ["soft / round", "hard / sharp"],
  },
];

function selectedImplement(value) {
  return implementFamilies.reduce((best, item) =>
    Math.abs(item.value - value) < Math.abs(best.value - value) ? item : best,
  );
}

export class PerformanceControls {
  constructor({ state, audition, scheduleRender, setStatus }) {
    this.state = state;
    this.audition = audition;
    this.scheduleRender = scheduleRender;
    this.setStatus = setStatus;
  }

  bind() {
    this.#bindMetallicControls();
    this.#bindKickControls();
    this.#bindMembraneControls();
  }

  paint() {
    if (this.state.recipeKey === "metal.cymbal.v1") this.#paintMetallic();
    else if (this.state.recipeKey === "drum.membrane.v1" ||
             this.state.recipeKey === "drum.snare.v1")
      this.#paintMembrane();
    else this.#paintKick();
  }

  #bindMetallicControls() {
    this.#bindSlider("hardness", "hardness", () =>
      this.state.eventDefaults.hardness);
    this.#bindSlider("contact-spread", "contactSpread", () =>
      this.state.eventDefaults.contactSpread);
    this.#bindSlider(
      "metal-constraint", "constraint",
      () => this.state.eventDefaults.constraint,
      amount => this.audition.setConstraint(amount),
    );
    document.querySelectorAll('input[name="implement"]').forEach(input => {
      input.onchange = event => {
        this.state.event.implement = Number(event.currentTarget.value);
        this.#paintMetallic();
        this.scheduleRender();
      };
    });
    byId("strike-pad").onpointerdown = event => {
      const bounds = event.currentTarget.getBoundingClientRect();
      this.state.event.location = clamp01(
        (event.clientX - bounds.left) / bounds.width);
      this.state.event.strength = Math.max(.02, clamp01(
        1 - (event.clientY - bounds.top) / bounds.height));
      ++this.state.event.seed;
      this.#trigger();
    };
  }

  #bindKickControls() {
    byId("kick-hardness").oninput = event => {
      this.state.event.hardness = Number(event.currentTarget.value);
      this.#paintKick();
      this.scheduleRender();
    };
    byId("kick-hardness").ondblclick = event => {
      event.preventDefault();
      event.currentTarget.value = .5;
      event.currentTarget.dispatchEvent(new Event("input"));
    };
    byId("kick-pad").onpointerdown = event => {
      const bounds = event.currentTarget.getBoundingClientRect();
      this.state.event.hardness = clamp01(
        (event.clientX - bounds.left) / bounds.width);
      this.state.event.strength = Math.max(.02, clamp01(
        1 - (event.clientY - bounds.top) / bounds.height));
      ++this.state.event.seed;
      this.#paintKick();
      this.#trigger();
    };
  }

  #bindMembraneControls() {
    this.#bindSlider("membrane-hardness", "hardness", () => .5);
    this.#bindSlider("membrane-contact-spread", "contactSpread", () => .25);
    document.querySelectorAll('input[name="membrane-implement"]').forEach(input => {
      input.onchange = event => {
        this.state.event.implement = Number(event.currentTarget.value);
        this.#paintMembrane();
        this.scheduleRender();
      };
    });
    byId("membrane-pad").onpointerdown = event => {
      const bounds = event.currentTarget.getBoundingClientRect();
      this.state.event.location = clamp01(
        (event.clientX - bounds.left) / bounds.width);
      this.state.event.strength = Math.max(.02, clamp01(
        1 - (event.clientY - bounds.top) / bounds.height));
      ++this.state.event.seed;
      this.#trigger();
    };
  }

  #bindSlider(id, key, defaultValue, onInput = () => {}) {
    byId(id).oninput = event => {
      this.state.event[key] = Number(event.currentTarget.value);
      event.currentTarget.nextElementSibling.textContent =
        this.state.event[key].toFixed(2);
      onInput(this.state.event[key]);
      this.scheduleRender();
    };
    byId(id).ondblclick = event => {
      event.preventDefault();
      event.currentTarget.value = defaultValue();
      event.currentTarget.dispatchEvent(new Event("input"));
    };
  }

  #paintMetallic() {
    const family = selectedImplement(this.state.event.implement);
    this.state.event.implement = family.value;
    document.querySelectorAll('input[name="implement"]').forEach(input => {
      input.checked = Number(input.value) === family.value;
    });
    byId("character-label").textContent = family.label;
    const endpoints = byId("character-endpoints").querySelectorAll("i");
    [endpoints[0].textContent, endpoints[1].textContent] = family.endpoints;
    this.#paintSlider("hardness", this.state.event.hardness);
    this.#paintSlider("contact-spread", this.state.event.contactSpread);
    this.#paintSlider("metal-constraint", this.state.event.constraint);
    this.audition.setConstraint(this.state.event.constraint);
  }

  #paintKick() {
    this.#paintSlider("kick-hardness", this.state.event.hardness);
  }

  #paintMembrane() {
    const family = selectedImplement(this.state.event.implement);
    this.state.event.implement = family.value;
    document.querySelectorAll('input[name="membrane-implement"]').forEach(input => {
      input.checked = Number(input.value) === family.value;
    });
    byId("membrane-character-label").textContent = family.label;
    const endpoints = byId("membrane-character-endpoints").querySelectorAll("i");
    [endpoints[0].textContent, endpoints[1].textContent] = family.endpoints;
    this.#paintSlider("membrane-hardness", this.state.event.hardness);
    this.#paintSlider("membrane-contact-spread", this.state.event.contactSpread);
  }

  #paintSlider(id, value) {
    byId(id).value = value;
    byId(id).nextElementSibling.textContent = value.toFixed(2);
  }

  #trigger() {
    this.audition.trigger({ ...this.state.event }).catch(
      error => this.setStatus(String(error)));
  }
}
