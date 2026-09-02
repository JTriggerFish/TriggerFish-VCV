import { copyFile, mkdir } from "node:fs/promises";
import { basename, dirname, join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const source = dirname(fileURLToPath(import.meta.url));
const entries = [];
for (let index = 2; index < process.argv.length; index += 2) {
  entries.push([process.argv[index].slice(2), process.argv[index + 1]]);
}
const options = Object.fromEntries(entries);
if (!options.module || !options.output) {
  throw new Error("expected --module and --output");
}

const output = resolve(options.output);
await mkdir(output, { recursive: true });
for (const name of [
  "index.html", "styles.css", "app.mjs", "fit_controls.mjs", "point_editor.mjs",
  "size_meta.mjs", "engine.mjs", "audio.mjs",
  "references.mjs", "analysis.mjs", "analysis_worker.mjs",
  "render_worker.mjs",
  "spectrogram.mjs", "state.mjs", "lookahead_limiter_processor.mjs",
  "reference_browser.mjs", "realtime_engine_worker.mjs",
  "ring_player_processor.mjs", "level_match.mjs",
]) {
  await copyFile(join(source, name), join(output, name));
}
const modulePath = resolve(options.module);
await copyFile(modulePath, join(output, basename(modulePath)));
await copyFile(
  modulePath.replace(/\.mjs$/, ".wasm"),
  join(output, basename(modulePath).replace(/\.mjs$/, ".wasm")),
);
await mkdir(join(output, "vendor"), { recursive: true });
await copyFile(
  join(source, "node_modules", "plotly.js-dist-min", "plotly.min.js"),
  join(output, "vendor", "plotly.min.js"),
);
console.log(`assembled ${output}`);
