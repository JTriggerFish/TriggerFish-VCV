import assert from "node:assert/strict";
import { svgPosition, decayDragDelta, shiftDecayPoints } from "../web/decay_curve_geometry.mjs";

// 600×220 content letterboxed inside a 200×220 CSS rectangle.
const point = svgPosition({ a: 1/3, b: 0, c: 0, d: 1/3, e: 10, f: 80 }, 110, 90);
assert.ok(Math.abs(point.x - 300) < 1e-10);
assert.ok(Math.abs(point.y - 30) < 1e-10);
assert.equal(svgPosition(null, 0, 0), null);
assert.equal(svgPosition({ a: 0, b: 0, c: 0, d: 0, e: 0, f: 0 }, 0, 0), null);
const normal = decayDragDelta(100, 90, -6, 5, 168, false);
const fine = decayDragDelta(100, 90, -6, 5, 168, true);
assert.equal(fine, normal*.1);
assert.ok(decayDragDelta(100, 100, -6, 5, 168, false) === 0);
assert.equal(normal + decayDragDelta(90, 100, -6, 5, 168, false), 0);
const points = [{ x: 40, y: 4 }, { x: 15000, y: 1 }];
const shifted = shiftDecayPoints(points, 3, -6, 5);
assert.deepEqual(shifted, [{ x: 40, y: 5 }, { x: 15000, y: 2 }]);
assert.equal(shifted[0].y - shifted[1].y, points[0].y - points[1].y);
assert.equal(points[0].y, 4);
console.log("Decay curve coordinate, precision and parallel-shift checks passed");
