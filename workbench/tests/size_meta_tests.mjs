import assert from "node:assert/strict";

import { expandedSizeMeta } from "../web/size_meta.mjs";

const descriptors = [
  { key: "impact_width", defaultValue: 1, minimum: .25, maximum: 4,
    scale: "logarithmic", index: 0 },
  { key: "bloom_energy_dependence", defaultValue: .35, minimum: 0, maximum: 1,
    scale: "linear", index: 1 },
  ...Array.from({ length: 6 }, (_, index) => ({
    key: `body_decay_active_${index + 1}`, defaultValue: 0,
    minimum: 0, maximum: 1, scale: "boolean", index: index + 2,
  })),
  ...Array.from({ length: 24 }, (_, index) => ({
    key: `resolved_frequency_${index}`, defaultValue: 100 * (index + 1),
    minimum: 40, maximum: 22000, scale: "logarithmic", index: index + 8,
  })),
  ...Array.from({ length: 24 }, (_, index) => ({
    key: `resolved_level_${index}`, defaultValue: 0,
    minimum: -24, maximum: 24, scale: "linear", index: index + 32,
  })),
];

const values = position => Object.fromEntries(
  expandedSizeMeta(descriptors, position).map(
    item => [item.descriptor.key, item.value],
  ),
);

assert.equal(values(.5).impact_width, 1);
assert.equal(values(.5).bloom_energy_dependence, .35);
assert.equal(values(0).impact_width, 1.8);
assert.equal(values(1).impact_width, .65);
assert.equal(values(0).resolved_frequency_5, 600 * .42);
assert.equal(values(1).resolved_frequency_5, 600 * 1.45);
assert.equal(values(0).resolved_frequency_23, 2400 * .42);
assert.equal(values(1).resolved_frequency_23, 2400 * 1.45);
assert.equal(values(0).resolved_level_0, 5);
assert.equal(values(1).resolved_level_23, -2);
for (let index = 1; index <= 6; ++index) {
  assert.equal(values(0)[`body_decay_active_${index}`], 0);
  assert.equal(values(1)[`body_decay_active_${index}`], 0);
}
assert.equal(values(0).dense_wash_frequency_3, undefined);
assert.equal(values(1).dense_wash_level_3, undefined);
assert.ok(Math.abs(values(.25).impact_width - Math.sqrt(1.8)) < 1.e-12);

console.log("size meta tests passed");
