import { MidiInputs } from "./midi_input.mjs";

const StorageKey = "tf-workbench-settings-v1";

function loadSettings() {
  try { return JSON.parse(localStorage.getItem(StorageKey)) ?? {}; }
  catch { return {}; }
}

export class SettingsController {
  constructor({ audition, onMidiNote, onStatus = () => {} }) {
    this.audition = audition;
    this.onStatus = onStatus;
    this.saved = loadSettings();
    this.midi = new MidiInputs({
      onNote: onMidiNote,
      onDevices: inputs => this.paintMidiInputs(inputs),
    });
  }

  bind() {
    this.dialog = document.getElementById("settings-dialog");
    this.input = document.getElementById("midi-input");
    this.channel = document.getElementById("midi-channel");
    document.getElementById("settings-open").onclick = () =>
      this.dialog.showModal();
    document.getElementById("settings-close").onclick = () => this.dialog.close();
    this.dialog.onclick = event => {
      if (event.target === this.dialog) this.dialog.close();
    };

    this.channel.value = String(this.saved.midiChannel ?? 0);
    this.midi.setChannel(this.channel.value);
    this.channel.onchange = () => {
      this.midi.setChannel(this.channel.value);
      this.save();
    };
    this.input.onchange = () => {
      this.midi.setInput(this.input.value);
      this.save();
    };
    const enable = document.getElementById("midi-enable");
    enable.disabled = !this.midi.supported;
    if (!this.midi.supported) {
      document.getElementById("midi-status").textContent =
        "Web MIDI is unavailable in this browser.";
    }
    enable.onclick = () => this.enableMidi(enable);
    if (this.midi.supported) {
      // MIDI device discovery is useful by default and does not need to start
      // the audio graph. Browsers that require an explicit permission gesture
      // leave the button available so the same request can be retried there.
      void this.enableMidi(enable, false);
    }
    this.paintAudioStatus();
  }

  async enableMidi(button, activateAudio = true) {
    button.disabled = true;
    document.getElementById("midi-status").textContent = "Requesting access…";
    try {
      await this.midi.enable();
      if (activateAudio) await this.audition.activate();
      button.textContent = "MIDI enabled";
    } catch (error) {
      button.disabled = false;
      document.getElementById("midi-status").textContent = String(error);
      this.onStatus(String(error));
    }
  }

  paintMidiInputs(inputs) {
    const selected = this.saved.midiInput ?? "all";
    this.input.replaceChildren(new Option("All connected inputs", "all"));
    for (const input of inputs) {
      const label = [input.manufacturer, input.name].filter(Boolean).join(" · ") ||
        "Unnamed MIDI input";
      this.input.add(new Option(label, input.id));
    }
    this.input.value = [...this.input.options].some(option => option.value === selected)
      ? selected : "all";
    this.input.disabled = false;
    this.channel.disabled = false;
    this.midi.setInput(this.input.value);
    document.getElementById("midi-status").textContent = inputs.length
      ? `${inputs.length} MIDI input${inputs.length === 1 ? "" : "s"} available`
      : "MIDI enabled; no input devices are connected.";
  }

  paintAudioStatus() {
    const latency = this.audition.latencyBreakdown;
    const rate = this.audition.sampleRate || 48000;
    const quantum = latency.workletQuantumMs || 1000 * 128 / rate;
    const parts = [
      `${latency.workletQuantumFrames}-frame direct worklet (${quantum.toFixed(1)} ms)`,
      `${latency.limiterMs.toFixed(1)} ms limiter`];
    if (latency.browserDeviceMs) {
      parts.push(`${latency.browserDeviceMs.toFixed(1)} ms browser/device`,
        `approximately ${latency.totalMs.toFixed(0)} ms total`);
    }
    document.getElementById("audio-buffer-status").textContent = parts.join(" · ");
  }

  save() {
    const value = {
      midiInput: this.input.value,
      midiChannel: Number(this.channel.value),
    };
    localStorage.setItem(StorageKey, JSON.stringify(value));
    this.saved = value;
  }
}
