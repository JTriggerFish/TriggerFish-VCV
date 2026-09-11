import assert from "node:assert/strict";
import { svgPosition, decayDragDelta, shiftDecayPoints, decayPosition, decaySeconds } from "../web/decay_curve_geometry.mjs";
import { DecayCurveEditor } from "../web/decay_curve_editor.mjs";

// 600×220 content letterboxed inside a 200×220 CSS rectangle.
const point = svgPosition({ a: 1/3, b: 0, c: 0, d: 1/3, e: 10, f: 80 }, 110, 90);
assert.ok(Math.abs(point.x - 300) < 1e-10);
assert.ok(Math.abs(point.y - 30) < 1e-10);
assert.equal(svgPosition(null, 0, 0), null);
assert.equal(svgPosition({ a: 0, b: 0, c: 0, d: 0, e: 0, f: 0 }, 0, 0), null);
const close = (a, b) => assert.ok(Math.abs(a-b)<1e-10, `${a} != ${b}`);
const normal = decayDragDelta(100, 90, -6, 5, 168, false, 3);
const fine = decayDragDelta(100, 90, -6, 5, 168, true, 3);
const position = value => decayPosition(2 ** value, 2 ** -6, 32);
close(position(3+fine)-position(3), (position(3+normal)-position(3))*.1);
close(decayDragDelta(100, 100, -6, 5, 168, false, 3), 0);
close(normal + decayDragDelta(90, 100, -6, 5, 168, false, 3+normal), 0);
for (const seconds of [.02, .1, 1, 3, 10, 20, 30])
  close(decaySeconds(decayPosition(seconds,.02,30),.02,30), seconds);
assert.ok(1-decayPosition(1,.02,30)>.79, 'Long decays get most of the display');
close(decayDragDelta(100, -1000, -6, 5, 168, false, 3), 2);
close(decayDragDelta(100, 1000, -6, 5, 168, false, 3), -9);
const points = [{ x: 40, y: 4 }, { x: 15000, y: 1 }];
const shifted = shiftDecayPoints(points, 3, -6, 5);
assert.deepEqual(shifted, [{ x: 40, y: 5 }, { x: 15000, y: 2 }]);
assert.equal(shifted[0].y - shifted[1].y, points[0].y - points[1].y);
assert.equal(points[0].y, 4);
// Displaying a softened log axis must not turn DSP interpolation into a straight line.
const editor = Object.create(DecayCurveEditor.prototype);
editor.width = 320;
editor.options = {minimumFrequency:40, maximumFrequency:15000,
  minimumLogSeconds:-6, maximumLogSeconds:5};
const coordinates = editor.curveCoordinates(points).split(' ').map(pair=>pair.split(',').map(Number));
assert.ok(coordinates.length>50);
for (const [x,y] of coordinates) {
  const fraction=(x-editor.xPosition(40))/(editor.xPosition(15000)-editor.xPosition(40));
  close(editor.logSeconds(y), 4-3*fraction);
}
close(editor.middleLevel(points), 2.5);
console.log("Decay curve coordinate, precision and parallel-shift checks passed");
