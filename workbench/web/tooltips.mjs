export class Tooltips {
  constructor() {
    this.popup = document.createElement("div");
    this.popup.id = "workbench-tooltip";
    this.popup.className = "tooltip";
    this.popup.setAttribute("role", "tooltip");
    document.body.append(this.popup);
    document.addEventListener("pointerover", event => this.enter(event));
    document.addEventListener("pointerout", event => this.leave(event));
    document.addEventListener("focusin", event => this.enter(event));
    document.addEventListener("focusout", event => this.leave(event));
    document.addEventListener("keydown", event => {
      if (event.key === "Escape") this.hide();
    });
  }

  target(event) {
    return event.target.closest?.("[data-tooltip]");
  }

  enter(event) {
    const target = this.target(event);
    if (!target || target === this.active) return;
    clearTimeout(this.timer);
    this.active = target;
    this.timer = setTimeout(() => this.show(target), 220);
  }

  leave(event) {
    const target = this.target(event);
    if (!target || target !== this.active || target.contains(event.relatedTarget)) {
      return;
    }
    this.hide();
  }

  show(target) {
    if (target !== this.active || !target.isConnected) return;
    this.popup.textContent = target.dataset.tooltip;
    target.setAttribute("aria-describedby", this.popup.id);
    const bounds = target.getBoundingClientRect();
    const width = Math.min(320, window.innerWidth - 20);
    this.popup.style.maxWidth = `${width}px`;
    this.popup.style.left = `${Math.max(10, Math.min(
      bounds.left, window.innerWidth - width - 10,
    ))}px`;
    this.popup.style.top = `${Math.min(
      window.innerHeight - 80, bounds.bottom + 7,
    )}px`;
    this.popup.classList.add("visible");
  }

  hide() {
    clearTimeout(this.timer);
    this.active?.removeAttribute("aria-describedby");
    this.active = null;
    this.popup.classList.remove("visible");
  }
}
