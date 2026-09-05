// Four-times, windowed-sinc reconstruction for the browser protection meter.
// The 32-frame analysis delay fits inside (does not extend) limiter lookahead.
export class TruePeakDetector {
  constructor() {
    this.radius = 32;
    this.size = 2 * this.radius + 1;
    this.ring = new Float64Array(this.size);
    this.write = 0;
    this.phases = [0.25, 0.5, 0.75].map(phase => {
      const coefficients = Float64Array.from({ length: this.size }, (_, tap) => {
        const x = tap - this.radius + phase;
        const sinc = Math.sin(Math.PI * x) / (Math.PI * x);
        const window = Math.abs(x) < this.radius
          ? 0.42 + 0.5 * Math.cos(Math.PI * x / this.radius) +
            0.08 * Math.cos(2 * Math.PI * x / this.radius) : 0;
        return sinc * window;
      });
      const sum = coefficients.reduce((a, b) => a + b, 0);
      return coefficients.map(value => value / sum);
    });
  }

  process(value) {
    this.ring[this.write] = value;
    let peak = Math.abs(value);
    for (const coefficients of this.phases) {
      let sum = 0;
      let index = this.write;
      for (const coefficient of coefficients) {
        sum += coefficient * this.ring[index];
        if (--index < 0) index = this.size - 1;
      }
      peak = Math.max(peak, Math.abs(sum));
    }
    if (++this.write === this.size) this.write = 0;
    return peak;
  }
}
