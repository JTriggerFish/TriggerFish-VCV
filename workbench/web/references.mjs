function text(view, offset, length) {
  return String.fromCharCode(
    ...new Uint8Array(view.buffer, view.byteOffset + offset, length),
  );
}

function pcmSample(view, offset, bits, format) {
  if (format === 3 && bits === 32) return view.getFloat32(offset, true);
  if (format !== 1) throw new Error(`unsupported WAV format ${format}`);
  if (bits === 16) return view.getInt16(offset, true) / 32768;
  if (bits === 24) {
    let value = view.getUint8(offset) | view.getUint8(offset + 1) << 8 |
      view.getUint8(offset + 2) << 16;
    if (value & 0x800000) value |= 0xff000000;
    return value / 8388608;
  }
  if (bits === 32) return view.getInt32(offset, true) / 2147483648;
  throw new Error(`unsupported ${bits}-bit WAV`);
}

export async function decodeWav(buffer, name, metadata = {}) {
  const digest = await crypto.subtle.digest("SHA-256", buffer);
  const sha256 = [...new Uint8Array(digest)]
    .map(value => value.toString(16).padStart(2, "0")).join("");
  const view = new DataView(buffer);
  if (text(view, 0, 4) !== "RIFF" || text(view, 8, 4) !== "WAVE") {
    throw new Error(`${name} is not a RIFF/WAVE file`);
  }
  let format;
  let dataOffset;
  let dataBytes;
  for (let offset = 12; offset + 8 <= view.byteLength;) {
    const kind = text(view, offset, 4);
    const size = view.getUint32(offset + 4, true);
    const payload = offset + 8;
    if (kind === "fmt ") {
      format = {
        tag: view.getUint16(payload, true),
        channels: view.getUint16(payload + 2, true),
        sampleRate: view.getUint32(payload + 4, true),
        blockAlign: view.getUint16(payload + 12, true),
        bits: view.getUint16(payload + 14, true),
      };
      if (format.tag === 0xfffe && size >= 40) {
        format.tag = view.getUint16(payload + 24, true);
      }
    } else if (kind === "data") {
      dataOffset = payload;
      dataBytes = size;
    }
    offset = payload + size + (size & 1);
  }
  if (!format || dataOffset === undefined) throw new Error("incomplete WAV file");
  const frames = Math.floor(dataBytes / format.blockAlign);
  const bytes = format.bits / 8;
  const samples = new Float32Array(frames);
  for (let frame = 0; frame < frames; ++frame) {
    let sum = 0;
    for (let channel = 0; channel < format.channels; ++channel) {
      sum += pcmSample(
        view, dataOffset + frame * format.blockAlign + channel * bytes,
        format.bits, format.tag,
      );
    }
    samples[frame] = sum / format.channels;
  }
  return {
    id: `sha256:${sha256}`,
    sha256,
    name,
    samples,
    sampleRate: format.sampleRate,
    channels: format.channels,
    bits: format.bits,
    duration: frames / format.sampleRate,
    ...metadata,
  };
}

export async function readWav(file) {
  return decodeWav(await file.arrayBuffer(), file.name);
}

export async function readRemoteReference(corpus, cell,
    referenceGainDb = corpus.reference_gain_db ?? 0) {
  const response = await fetch(cell.url);
  if (!response.ok) throw new Error(`could not load ${cell.label}`);
  const reference = await decodeWav(await response.arrayBuffer(), cell.label, {
    corpus: {
      id: corpus.id, name: corpus.name,
      referenceGainDb: corpus.reference_gain_db ?? 0,
    },
    cell: { ...cell },
  });
  return calibrateReference(reference, referenceGainDb);
}

// One fixed corpus conversion, not per-hit normalization. Keep source hashes
// untouched and record the conversion so snapshots/analysis remain auditable.
export function calibrateReference(reference, gainDb) {
  if (!Number.isFinite(gainDb) || gainDb < -120 || gainDb > 120)
    throw new Error("invalid reference calibration gain");
  if (reference.referenceGainDb !== undefined)
    throw new Error("reference is already calibrated");
  return setReferenceGain(reference, gainDb);
}

// Always derive manual edits from immutable source samples, never a previous
// scaled buffer. This keeps reset and snapshot restore independent of edit order.
export function setReferenceGain(reference, gainDb) {
  if (!Number.isFinite(gainDb) || gainDb < -120 || gainDb > 120)
    throw new Error("invalid reference gain");
  const rawSamples = reference.rawSamples ?? reference.samples;
  const gain = 10 ** (gainDb / 20);
  return {
    ...reference, rawSamples, referenceGainDb: gainDb,
    samples: rawSamples.map(sample => sample * gain),
  };
}

export async function readReferenceCorpora() {
  const response = await fetch("/api/reference-corpora");
  if (!response.ok) throw new Error("could not read reference corpus index");
  return (await response.json()).corpora;
}

export async function readReferences(files) {
  const results = [];
  for (const file of files) results.push(await readWav(file));
  return results;
}
