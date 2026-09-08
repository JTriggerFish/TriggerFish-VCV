import assert from "node:assert/strict";
import { mountDiagnosticAudition } from "../web/diagnostic_audition.mjs";

// No audio context/device: verify routing through the supplied safe player.
let panel, played, status, requests = 0;
const element = () => ({ children: [], append(...items) { this.children.push(...items); } });
globalThis.document = {
  createElement: element,
  querySelector: () => ({ after(value) { panel = value; } }),
};
globalThis.location = new URL("http://localhost/");
const wav = new ArrayBuffer(48);
const view = new DataView(wav);
const text = (offset, value) => [...value].forEach((c, i) => view.setUint8(offset+i, c.charCodeAt(0)));
text(0, "RIFF"); view.setUint32(4, 40, true); text(8, "WAVE");
text(12, "fmt "); view.setUint32(16, 16, true);
view.setUint16(20, 3, true); view.setUint16(22, 1, true);
view.setUint32(24, 48000, true); view.setUint32(28, 192000, true);
view.setUint16(32, 4, true); view.setUint16(34, 32, true);
text(36, "data"); view.setUint32(40, 4, true); view.setFloat32(44, .25, true);
globalThis.fetch = async url => {
  requests++;
  return String(url).endsWith("manifest.json")
    ? { ok: true, json: async () => ({ title: "Rate check", note: "Fixed controls",
      clips: [{ label: "4x single", file: "test.wav" }] }) }
    : { ok: true, arrayBuffer: async () => wav };
};
const audition = { play: async (samples, rate) => { played = { samples, rate }; } };
await mountDiagnosticAudition(audition, value => { status = value; });
assert.equal(requests, 0, "default workbench does not load diagnostics");
globalThis.location = new URL("http://localhost/?audition=checks/manifest.json");
await mountDiagnosticAudition(audition, value => { status = value; });
assert.equal(panel.children[0].textContent, "Rate check");
assert.equal(played, undefined, "never autoplay");
await panel.children[1].onclick();
assert.equal(played.rate, 48000);
assert.equal(played.samples[0], .25, "no hidden gain/normalization");
assert.match(status, /4x single/);
globalThis.location = new URL("http://localhost/?audition=https://elsewhere.test/a.json");
await mountDiagnosticAudition(audition, value => { status = value; });
assert.ok(status instanceof Error, "diagnostic failures reach the persistent error banner");
assert.match(status.message, /same-origin/);
assert.equal(requests, 2);
console.log("Diagnostic audition tests passed");
