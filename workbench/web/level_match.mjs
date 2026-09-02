function meanSquare(samples, sampleRate, seconds) {
  const count = Math.min(
    samples?.length ?? 0,
    Math.max(1, Math.round(sampleRate * seconds)),
  );
  if (!count) return 0;
  let energy = 0;
  for (let index = 0; index < count; ++index) {
    const sample = Number(samples[index]);
    if (Number.isFinite(sample)) energy += sample * sample;
  }
  return energy / count;
}

export function comparisonLevelDb(samples, sampleRate, seconds = 1) {
  return 10 * Math.log10(Math.max(1e-20, meanSquare(
    samples, sampleRate, seconds,
  )));
}

export function matchedModelLevelDb({
  currentDb, reference, referenceSampleRate, synthesis, synthesisSampleRate,
  minimumDb, maximumDb, seconds = 1,
}) {
  const referenceDb = comparisonLevelDb(
    reference, referenceSampleRate, seconds,
  );
  const synthesisDb = comparisonLevelDb(
    synthesis, synthesisSampleRate, seconds,
  );
  return Math.max(minimumDb, Math.min(
    maximumDb, currentDb + referenceDb - synthesisDb,
  ));
}
