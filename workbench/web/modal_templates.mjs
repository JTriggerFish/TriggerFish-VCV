// Design tools only: these write ordinary editable modes, never runtime rules.
// Lowest 16 distinct roots j_(m,n) of J_m, divided by j_(0,1).
// Generated offline with scipy.special.jn_zeros (m=0..11, n=1..8), sorted.
// Angular degeneracy is not duplicated; these are editable design handles.
export const MembraneRatios = Object.freeze([
  1, 1.593340506, 2.135548787, 2.295417267,
  2.653066405, 2.917295455, 3.155464815, 3.500147490,
  3.598484674, 3.647451179, 4.058931883, 4.131738160,
  4.230439128, 4.601044534, 4.610051645, 4.831885263,
]);

// C1-continuous bend above a protected low core. Independent of total count:
// adding more modes must not retune existing centres. Not a physical gong law.
function stretchedRatio(family, index, stretch, harmonicCore) {
  const ratio = family === "membrane" ? MembraneRatios[index] : index + 1;
  const upper = Math.max(0, (index + 1 - harmonicCore) / harmonicCore);
  return ratio * Math.hypot(1, stretch * upper);
}

// Inverse used by offline fitting to request a top frequency using the SAME law.
export function modalTemplateStretch({family = "harmonic", fundamental, count,
  topFrequency, harmonicCore = 4}) {
  if (![fundamental, count, topFrequency, harmonicCore].every(Number.isFinite) ||
      !["harmonic", "membrane"].includes(family) || fundamental <= 0 ||
      !Number.isInteger(count) || count < 1 || count > (family === "membrane" ? 16 : 32) ||
      !Number.isInteger(harmonicCore) || harmonicCore < 1 || harmonicCore > 8)
    throw Error("Invalid stretch endpoint settings");
  const base = fundamental * stretchedRatio(family, count - 1, 0, harmonicCore);
  const upper = Math.max(0, (count - harmonicCore) / harmonicCore);
  if (Math.abs(topFrequency - base) < base * 1e-12) return 0;
  const stretch = Math.sqrt((topFrequency / base) ** 2 - 1) / upper;
  if (topFrequency < base || !Number.isFinite(stretch) || stretch > 1 + 1e-12)
    throw Error("Top frequency is outside this harmonic core/stretch range");
  return Math.min(1, stretch);
}

export function modalTemplate({
  family = "membrane", fundamental = 55, count = 16,
  level = 0, rolloff = 6, turbulence = 1, stretch = 0, harmonicCore = 4,
  minimumFrequency = 20, maximumFrequency = 15000,
} = {}) {
  if (!["membrane", "harmonic"].includes(family) ||
      ![fundamental, count, level, rolloff, turbulence, stretch, harmonicCore, minimumFrequency, maximumFrequency].every(Number.isFinite) ||
      fundamental <= 0 || minimumFrequency <= 0 || maximumFrequency < minimumFrequency ||
      turbulence < 0 || turbulence > 2 || stretch < 0 || stretch > 1 ||
      !Number.isInteger(harmonicCore) || harmonicCore < 1 || harmonicCore > 8 ||
      count < 1 || count > 32 || !Number.isInteger(count))
    throw Error("Invalid modal template settings");
  const limit = modalTemplateLimit({family, fundamental, stretch, harmonicCore, minimumFrequency, maximumFrequency});
  if (count > limit) throw Error(`Only ${limit} modes fit this formula and frequency range`);
  const length = count;
  return Array.from({length}, (_, index) => {
    const ratio = stretchedRatio(family, index, stretch, harmonicCore);
    const position = length > 1 ? index / (length - 1) : 0;
    return {
      frequency: fundamental * ratio,
      level: Math.max(-72, Math.min(6, level - rolloff * Math.log2(ratio))),
      centre: index === 0 ? 1 : .65 * (index % 2 ? -1 : 1) * (1 - .45 * position),
      edge: .18 + .82 * position,
      turbulence, active: level - rolloff * Math.log2(ratio) > -72,
    };
  });
}

// Limits are shared by the UI preview and the pure generator. Never truncate.
export function modalTemplateLimit({family, fundamental, stretch = 0, harmonicCore = 4, minimumFrequency = 20,
  maximumFrequency = 15000, capacity = 32}) {
  if (![fundamental, stretch, harmonicCore, minimumFrequency, maximumFrequency, capacity].every(Number.isFinite) ||
      !["harmonic", "membrane"].includes(family) || stretch < 0 || stretch > 1 ||
      !Number.isInteger(harmonicCore) || harmonicCore < 1 || harmonicCore > 8 ||
      minimumFrequency <= 0 || fundamental < minimumFrequency || fundamental > maximumFrequency ||
      !Number.isInteger(capacity) || capacity < 1 || capacity > 32) return 0;
  // Evaluate the same expression as generation, including boundary rounding.
  const ratios = family === "membrane" ? MembraneRatios : Array.from({length:32},(_,i)=>i+1);
  return Math.min(capacity, ratios.filter((_, index) =>
    fundamental * stretchedRatio(family, index, stretch, harmonicCore) <= maximumFrequency).length);
}
