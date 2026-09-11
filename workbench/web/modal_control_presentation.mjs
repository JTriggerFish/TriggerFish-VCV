// Presentation only: keys, values and DSP mappings stay unchanged.
export const shimmerPosition = (descriptor, value) => Math.sqrt(value / descriptor.maximum);
export const shimmerValue = (descriptor, position) => descriptor.maximum * position ** 2;

export const modalControlTitles = {
  body_tune: "Tuning",
  field_turbulence: "Surrounding rings at 1 kHz",
  field_turbulence_slope: "Bass / treble balance",
  field_packet_spread: "Spread",
  field_satellite_density: "Density",
  field_doublet_split: "Speed at 125 Hz",
  field_beat_depth: "Depth",
  field_beat_rate_tilt: "Treble speed scaling",
  field_wander_hz: "Amount",
  field_wander_rate: "Speed",
  field_motion_depth: "Amount",
  field_motion_rate: "Speed",
  field_motion_sharing: "Shimmer moves together",
  field_phase_bandwidth: "Amount",
  field_phase_tilt: "Bass / treble balance",
};

export function modalControlActivity(value) {
  const paired = [2, 3, 4].includes(Math.round(value("field_distribution")));
  const beating = paired && value("field_beat_depth") > 0;
  return {
    field_beat_depth: paired,
    field_doublet_split: beating,
    field_beat_rate_tilt: beating && value("field_doublet_split") > 0,
    field_wander_rate: value("field_wander_hz") > 0,
    field_motion_rate: value("field_motion_depth") > 0,
    field_motion_sharing: value("field_motion_depth") > 0,
    field_phase_tilt: value("field_phase_bandwidth") > 0,
  };
}
