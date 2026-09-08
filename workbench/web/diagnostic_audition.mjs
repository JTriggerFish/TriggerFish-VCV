import { decodeWav } from "./references.mjs";

// Opt-in offline comparisons use the ordinary master/limiter, not HTML audio.
// They never replace the reference, patch, or current synthesis.
export async function mountDiagnosticAudition(audition, setStatus) {
  const path = new URLSearchParams(location.search).get("audition");
  if (!path) return;
  try {
    const manifestUrl = localUrl(path, location.href);
    const response = await fetch(manifestUrl, { cache: "no-store" });
    if (!response.ok) throw new Error(`Diagnostic manifest: HTTP ${response.status}`);
    const manifest = await response.json();
    const clips = await Promise.all(manifest.clips.map(async clip => {
      const response = await fetch(localUrl(clip.file, manifestUrl));
      if (!response.ok) throw new Error(`Diagnostic audio: HTTP ${response.status}`);
      const audio = await decodeWav(await response.arrayBuffer(), clip.label);
      return { ...audio, label: clip.label };
    }));
    const panel = document.createElement("section");
    panel.className = "toolbar";
    const title = document.createElement("strong");
    title.textContent = manifest.title;
    panel.append(title);
    for (const clip of clips) {
      const button = document.createElement("button");
      button.textContent = clip.label;
      button.onclick = () => audition.play(clip.samples, clip.sampleRate)
        .then(() => setStatus(`Diagnostic playback: ${clip.label}`))
        .catch(error => setStatus(String(error)));
      panel.append(button);
    }
    const note = document.createElement("span");
    note.textContent = manifest.note;
    panel.append(note);
    document.querySelector("nav.toolbar").after(panel);
  } catch (error) { setStatus(new Error(`Diagnostic audition: ${error.message}`)); }
}

function localUrl(path, base) {
  const url = new URL(path, base);
  if (url.origin !== location.origin) throw new Error("Use same-origin diagnostic files");
  return url;
}
