import assert from "node:assert/strict";
import {BloomTimingKeys, bloomTimingValues} from "../web/bloom_timing_meta.mjs";
import {bloomRateNormalized, bloomRateDenormalized} from "../web/bloom_control_scaling.mjs";

const descriptors = [
  {key:"bloom_rate", minimum:0, maximum:16},
  {key:"body_brightness", minimum:-72, maximum:24},
  {key:"body_excitation_centre", minimum:1, maximum:15000},
];
const original = {bloom_rate:2, body_brightness:-18, body_excitation_centre:1000};
assert.deepEqual(bloomTimingValues(original,0,descriptors).values, original);
for (const position of [-1,-.5,.5,1]) {
  const {values} = bloomTimingValues(original,position,descriptors);
  assert.deepEqual(Object.keys(values), BloomTimingKeys);
  for (const key of BloomTimingKeys)
    assert.equal(Math.sign(values[key]-original[key]), -Math.sign(position));
}
assert.equal(bloomTimingValues({...original,bloom_rate:12},-1,descriptors).limited.length, 1);
assert.deepEqual(original, {bloom_rate:2,body_brightness:-18,body_excitation_centre:1000});
// Returning from a clipped extreme uses the original, never accumulated deltas.
assert.deepEqual(bloomTimingValues(original,0,descriptors).values, original);
assert.equal(bloomTimingValues({...original,bloom_rate:0},1,descriptors).values.bloom_rate,0);
for (const invalid of [NaN, Infinity, -2, 2])
  assert.throws(()=>bloomTimingValues(original,invalid,descriptors));
assert.throws(()=>bloomTimingValues({...original,bloom_rate:-1},0,descriptors));
const split = {...original,bloom_energy_acceleration:.2,bloom_energy_sensitivity:1.3};
const moved = {...split,...bloomTimingValues(split,1,descriptors).values};
assert.equal(moved.bloom_energy_acceleration,.2);
assert.equal(moved.bloom_energy_sensitivity,1.3);
for(const value of [0,.01,.1,1,4,16]) {
  const position = bloomRateNormalized(descriptors[0],value);
  assert.ok(Math.abs(bloomRateDenormalized(descriptors[0],position)-value)<1e-10);
}
console.log("Bloom timing: three visible parameters; identity, bounds and reversibility pass");
