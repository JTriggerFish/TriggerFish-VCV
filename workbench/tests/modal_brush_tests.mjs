import assert from 'node:assert/strict';
import {ModalEditor} from '../web/modal_editor.mjs';

const erb = f => 21.4 * Math.log10(1 + .00437 * f);
const inverseErb = e => (10 ** (e / 21.4) - 1) / .00437;
const centre = erb(1000);
const editor = Object.create(ModalEditor.prototype);
editor.options = {minimumLevel: -72, maximumLevel: 6};
const points = levels => levels.map((level, i) => ({
  frequency: inverseErb(centre + i - 1), level, active: true,
}));

// A quieter centre must not pull down its louder neighbours when raised.
let bars = points([-10, -30, -10]);
editor.brushLevels(bars, centre, 1, -20);
assert.ok(bars[0].level > -10 && bars[1].level > -30 && bars[2].level > -10);
assert.ok(bars[1].level + 30 > bars[0].level + 10, 'Centre gets the strongest change');
assert.ok(Math.abs(bars[0].level - bars[2].level) < 1e-10, 'Falloff is symmetric');

// The reverse gesture must not raise quiet neighbours while lowering a peak.
bars = points([-40, -10, -40]);
editor.brushLevels(bars, centre, 1, -20);
assert.ok(bars[0].level < -40 && bars[1].level < -10 && bars[2].level < -40);

// Steady height must not keep flattening the surrounding shape.
bars = points([-10, -30, -10]);
const original = structuredClone(bars);
editor.brushLevels(bars, centre, 1, -30);
assert.deepEqual(bars, original);

// Locality, inactive handles and UI bounds.
bars = points([5, -30, -10]);
bars.push({frequency: 15000, level: -24, active: true},
  {frequency: 1000, level: -72, active: false});
editor.brushLevels(bars, centre, 1, 6);
assert.equal(bars[0].level, 6);
assert.equal(bars[3].level, -24);
assert.equal(bars[4].level, -72);
assert.ok(bars.every(p => p.level >= -72 && p.level <= 6));

// Both tools use the same operation; insertion keeps the clicked height.
for (const kind of ['shape', 'paint']) {
  editor.drag = {kind};
  editor.brushErb = 1;
  editor.frequency = () => 1000;
  editor.level = () => -20;
  editor.snapFrequency = f => f;
  bars = points([-10, -30, -10]);
  editor.paintAt({x: 0, y: 0}, {}, bars);
  assert.ok(bars[0].level > -10 && bars[2].level > -10);
}
bars = [{active: false}];
editor.paintAt({x: 0, y: 0}, {}, bars);
assert.equal(bars[0].level, -20);
assert.equal(bars[0].frequency, 1000);
console.log('Modal brush: direction, shape, locality, bounds and insertion pass');
