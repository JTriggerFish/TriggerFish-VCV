import { copyFile, mkdir } from "node:fs/promises";
import { basename, dirname, join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const source = dirname(fileURLToPath(import.meta.url));
const entries = [];
for (let index = 2; index < process.argv.length; index += 2) {
  entries.push([process.argv[index].slice(2), process.argv[index + 1]]);
}
const options = Object.fromEntries(entries);
if (!options.module || !options["audio-module"] || !options.output) {
  throw new Error("expected --module, --audio-module and --output");
}

const output = resolve(options.output);
await mkdir(output, { recursive: true });
for (const name of [
  "index.html", "styles.css", "app.mjs", "fit_controls.mjs", "point_editor.mjs",
  "percussion_patch.mjs", "percussion_registry.mjs",
  "metallic_plate_patch.mjs", "compact_kick_patch.mjs", "membrane_patch.mjs",
  "recipe_adapter.mjs", "routing_view.mjs", "kick_controls.mjs",
  "membrane_controls.mjs",
  "routing_controller.mjs", "performance_controls.mjs",
  "analysis_controls.mjs", "recipe_controller.mjs", "waveform_view.mjs",
  "modal_editor.mjs", "decay_curve_editor.mjs", "tooltips.mjs",
  "size_meta.mjs", "engine.mjs", "wasm_engine_core.mjs", "audio.mjs",
  "standby_renderer.mjs", "configuration_preparer.mjs",
  "limiter_config.mjs", "midi_input.mjs", "settings.mjs",
  "references.mjs", "analysis.mjs", "analysis_worker.mjs",
  "render_worker.mjs", "preparation_worker.mjs",
  "spectrogram.mjs", "state.mjs", "lookahead_limiter_processor.mjs",
  "reference_browser.mjs", "percussion_audio_worklet_processor.mjs",
  "level_match.mjs",
]) {
  await copyFile(join(source, name), join(output, name));
}
const modulePath = resolve(options.module);
await copyFile(modulePath, join(output, basename(modulePath)));
await copyFile(
  modulePath.replace(/\.mjs$/, ".wasm"),
  join(output, basename(modulePath).replace(/\.mjs$/, ".wasm")),
);
const audioModulePath = resolve(options["audio-module"]);
await copyFile(audioModulePath, join(output, basename(audioModulePath)));
await mkdir(join(output, "vendor"), { recursive: true });
await copyFile(
  join(source, "node_modules", "plotly.js-dist-min", "plotly.min.js"),
  join(output, "vendor", "plotly.min.js"),
);
console.log(`assembled ${output}`);
