function meanSquare(samples, sampleRate, seconds, startSeconds = 0) {
  const first = Math.min(
    samples?.length ?? 0, Math.max(0, Math.round(sampleRate * startSeconds)),
  );
  const count = Math.min(
    (samples?.length ?? 0) - first,
    Math.max(1, Math.round(sampleRate * seconds)),
  );
  if (!count) return 0;
  let energy = 0;
  for (let index = 0; index < count; ++index) {
    const sample = Number(samples[first + index]);
    if (Number.isFinite(sample)) energy += sample * sample;
  }
  return energy / count;
}

export function comparisonLevelDb(
  samples, sampleRate, seconds = .3, startSeconds = 0,
) {
  return 10 * Math.log10(Math.max(1e-20, meanSquare(
    samples, sampleRate, seconds, startSeconds,
  )));
}

export function matchedModelLevelDb({
  currentDb, reference, referenceSampleRate, synthesis, synthesisSampleRate,
  minimumDb, maximumDb, seconds = .3,
  referenceStartSeconds = 0, synthesisStartSeconds = 0,
}) {
  const referenceDb = comparisonLevelDb(
    reference, referenceSampleRate, seconds, referenceStartSeconds,
  );
  const synthesisDb = comparisonLevelDb(
    synthesis, synthesisSampleRate, seconds, synthesisStartSeconds,
  );
  return Math.max(minimumDb, Math.min(
    maximumDb, currentDb + referenceDb - synthesisDb,
  ));
}
