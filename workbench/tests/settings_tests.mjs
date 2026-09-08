import assert from "node:assert/strict";

import { decodeNoteOn, MidiInputs } from "../web/midi_input.mjs";

assert.deepEqual(decodeNoteOn(new Uint8Array([0x90, 60, 127])), {
  note: 60, velocity: 1, channel: 1,
});
assert.deepEqual(decodeNoteOn(new Uint8Array([0x9f, 36, 64]), 16), {
  note: 36, velocity: 64 / 127, channel: 16,
});
assert.equal(decodeNoteOn(new Uint8Array([0x90, 60, 0])), null);
assert.equal(decodeNoteOn(new Uint8Array([0x80, 60, 100])), null);
assert.equal(decodeNoteOn(new Uint8Array([0x91, 60, 100]), 1), null);

const received = [];
let enumerated = [];
const first = { id: "first", state: "connected", name: "First" };
const second = { id: "second", state: "connected", name: "Second" };
const disconnected = { id: "gone", state: "disconnected", name: "Gone" };
const midi = new MidiInputs({
  onNote: event => received.push(event),
  onDevices: inputs => { enumerated = inputs; },
});
midi.access = { inputs: new Map([
  [first.id, first], [second.id, second], [disconnected.id, disconnected],
]) };
midi.refresh();
assert.deepEqual(enumerated, [first, second]);
first.onmidimessage({ data: new Uint8Array([0x90, 60, 100]) });
assert.equal(received.length, 1);
assert.equal(received[0].input, first);
midi.setInput("second");
first.onmidimessage({ data: new Uint8Array([0x90, 61, 100]) });
second.onmidimessage({ data: new Uint8Array([0x90, 62, 100]) });
assert.equal(received.length, 2);
assert.equal(received[1].note, 62);
midi.setChannel(2);
second.onmidimessage({ data: new Uint8Array([0x90, 63, 100]) });
second.onmidimessage({ data: new Uint8Array([0x91, 64, 100]) });
assert.equal(received.length, 3);
assert.equal(received[2].channel, 2);

// An expected denial of automatic discovery is a settings notice. An explicit
// user request and unrelated failures must still reach the global error UI.
const {SettingsController} = await import("../web/settings.mjs");
const status = {textContent:""};
globalThis.document = {getElementById:()=>status};
const reports = [];
const denied = new DOMException("Permission denied", "NotAllowedError");
const controller = {
  midi:{enable:async()=>{throw denied;}},
  audition:{activate:async()=>{}}, onStatus:error=>reports.push(error),
};
const button={disabled:false};
await SettingsController.prototype.enableMidi.call(controller,button,false);
assert.equal(reports.length,0);
assert.match(status.textContent,/site permissions/);
assert.equal(button.disabled,false);
await SettingsController.prototype.enableMidi.call(controller,button,true);
assert.equal(reports.length,1);
controller.midi.enable=async()=>{throw Error("Broken MIDI service");};
await SettingsController.prototype.enableMidi.call(controller,button,false);
assert.equal(reports.length,2);
controller.midi.enable=async()=>{};
controller.audition.activate=async()=>{throw denied;};
await SettingsController.prototype.enableMidi.call(controller,button,true);
assert.equal(reports.length,3); // Audio permission failure is not a MIDI notice.
delete globalThis.document;
console.log("workbench settings tests passed");
