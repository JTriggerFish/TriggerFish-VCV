import assert from "node:assert/strict";

import { expandedSizeMeta } from "../web/size_meta.mjs";

const descriptors = [
  { key: "impact_width", defaultValue: 1, minimum: .25, maximum: 4,
    scale: "logarithmic", index: 0 },
  { key: "bloom_nonlinearity", defaultValue: .35, minimum: 0, maximum: 1,
    scale: "linear", index: 1 },
  ...Array.from({ length: 12 }, (_, index) => ({
    key: `resolved_frequency_${index}`, defaultValue: 100 * (index + 1),
    minimum: 40, maximum: 22000, scale: "logarithmic", index: index + 2,
  })),
  ...Array.from({ length: 12 }, (_, index) => ({
    key: `resolved_level_${index}`, defaultValue: 0,
    minimum: -24, maximum: 24, scale: "linear", index: index + 14,
  })),
  ...Array.from({ length: 8 }, (_, index) => ({
    key: `dense_wash_frequency_${index}`, defaultValue: 200 * (index + 1),
    minimum: 40, maximum: 22000, scale: "logarithmic", index: index + 26,
  })),
  ...Array.from({ length: 8 }, (_, index) => ({
    key: `dense_wash_level_${index}`, defaultValue: 0,
    minimum: -24, maximum: 24, scale: "linear", index: index + 34,
  })),
];

const values = position => Object.fromEntries(
  expandedSizeMeta(descriptors, position).map(
    item => [item.descriptor.key, item.value],
  ),
);

assert.equal(values(.5).impact_width, 1);
assert.equal(values(.5).bloom_nonlinearity, .35);
assert.equal(values(0).impact_width, 1.8);
assert.equal(values(1).impact_width, .65);
assert.equal(values(0).resolved_frequency_5, 600 * .42);
assert.equal(values(1).resolved_frequency_5, 600 * 1.45);
assert.equal(values(0).resolved_level_0, 5);
assert.equal(values(1).resolved_level_5, 7);
assert.equal(values(0).dense_wash_frequency_3, 800 * .42);
assert.equal(values(1).dense_wash_level_3, 5);
assert.ok(Math.abs(values(.25).impact_width - Math.sqrt(1.8)) < 1.e-12);

console.log("size meta tests passed");
