export function decodeNoteOn(data, channelFilter = 0) {
  if (!data || data.length < 3) return null;
  const status = data[0];
  const channel = (status & 0x0f) + 1;
  if ((status & 0xf0) !== 0x90 || data[2] === 0) return null;
  if (channelFilter && channel !== channelFilter) return null;
  return { note: data[1], velocity: data[2] / 127, channel };
}

export class MidiInputs {
  constructor({ onNote = () => {}, onDevices = () => {} } = {}) {
    this.onNote = onNote;
    this.onDevices = onDevices;
    this.access = null;
    this.inputId = "all";
    this.channel = 0;
  }

  get supported() {
    return typeof navigator !== "undefined" &&
      typeof navigator.requestMIDIAccess === "function";
  }

  async enable() {
    if (!this.supported) throw new Error("Web MIDI is unavailable in this browser");
    this.access = await navigator.requestMIDIAccess({ sysex: false });
    this.access.onstatechange = () => this.refresh();
    this.refresh();
  }

  setInput(inputId) {
    this.inputId = inputId || "all";
  }

  setChannel(channel) {
    this.channel = Math.max(0, Math.min(16, Number(channel) || 0));
  }

  refresh() {
    const inputs = [...(this.access?.inputs.values() ?? [])]
      .filter(input => input.state === "connected");
    for (const input of inputs) {
      input.onmidimessage = event => {
        if (this.inputId !== "all" && input.id !== this.inputId) return;
        const note = decodeNoteOn(event.data, this.channel);
        if (note) this.onNote({ ...note, input });
      };
    }
    this.onDevices(inputs);
  }
}
