import assert from "node:assert/strict";
import {modalControlActivity, modalControlTitles, shimmerPosition, shimmerValue} from "../web/modal_control_presentation.mjs";

const values = {field_distribution:3, field_beat_depth:.3, field_doublet_split:1,
  field_wander_hz:0, field_motion_depth:0, field_phase_bandwidth:0};
const before = {...values};
const activity = () => modalControlActivity(key => values[key]);
assert.equal(activity().field_beat_depth, true);
assert.equal(activity().field_doublet_split, true);
for (const key of ["field_wander_rate", "field_motion_rate", "field_motion_sharing", "field_phase_tilt"])
  assert.equal(activity()[key], false);
assert.deepEqual(values, before, "presentation must not change DSP values");
values.field_wander_hz = .3;
values.field_motion_depth = 1.5;
values.field_phase_bandwidth = .2;
assert.ok(Object.values(activity()).every(Boolean));
values.field_beat_depth = 0;
assert.equal(activity().field_beat_depth, true, "depth can enable beating again");
assert.equal(activity().field_doublet_split, false);
assert.equal(activity().field_beat_rate_tilt, false);
for (const layout of [0, 1]) {
  values.field_distribution = layout;
  assert.equal(activity().field_beat_depth, false);
}
for (const layout of [2, 3, 4]) {
  values.field_distribution = layout;
  assert.equal(activity().field_beat_depth, true);
}
assert.equal(modalControlTitles.field_motion_depth, "Amount");
assert.equal(modalControlTitles.field_wander_hz, "Amount");
const descriptor = {maximum:12};
for (const value of [0, .001, .1, .6, 1.5, 6, 12])
  assert.ok(Math.abs(shimmerValue(descriptor, shimmerPosition(descriptor, value))-value)<1e-12);
assert.ok(Math.abs(shimmerValue(descriptor,.1)-.12)<1e-12);
assert.equal(shimmerValue(descriptor,1), 12);
console.log("Modal presentation preserves values, gates inactive dependents and refines shimmer travel");
