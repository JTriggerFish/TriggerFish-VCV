// Meter ballistics only: never changes the limiter or audio signal.
export class LimiterMeter {
  constructor() { this.reset(); }

  reset() {
    this.maximumDb = 0;
    this.recentDb = 0;
    this.activeUntil = 0;
  }

  ingest(message, now = performance.now()) {
    const db = Math.max(0, -(message.intervalReductionDb ?? message.reductionDb));
    if (!Number.isFinite(db)) return;
    this.maximumDb = Math.max(this.maximumDb, db);
    if (db >= .1) {
      this.recentDb = now < this.activeUntil ? Math.max(this.recentDb, db) : db;
      this.activeUntil = now + 1000;
    }
  }

  read(now = performance.now()) {
    return { active: now < this.activeUntil,
      recentDb: now < this.activeUntil ? this.recentDb : 0,
      maximumDb: this.maximumDb };
  }
}

export function paintLimiterMeter(element, reading) {
  element.classList.toggle('is-limiting', reading.active);
  element.classList.toggle('has-limited', reading.maximumDb >= .1);
  element.querySelector('.limiter-label').textContent = reading.active ? 'LIMITING' : 'Limiter';
  const bar = element.querySelector('meter');
  bar.value = Math.min(bar.max, reading.recentDb);
  bar.setAttribute('aria-valuetext', `${reading.recentDb.toFixed(1)} dB gain reduction`);
  element.querySelector('output').textContent = `${reading.recentDb.toFixed(1)} dB`;
  element.querySelector('button').textContent = `Max ${reading.maximumDb.toFixed(1)} dB · reset`;
}
