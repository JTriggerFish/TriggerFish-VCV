import assert from 'node:assert/strict';
import {packetAllocation} from '../web/packet_allocation.mjs';
import {ModalEditor} from '../web/modal_editor.mjs';
import {hasPairedRing,beatRatePosition,beatRateValue,ringBeatRate,beatDepthPosition,beatDepthValue} from '../web/packet_layout_control.mjs';

for (const layout of [3,4]) {
  const editor={options:{packetLayout:()=>layout}}, point={frequency:8000};
  const frequency=ModalEditor.prototype.packetFrequency.call(editor,point,4);
  const distance=ModalEditor.prototype.packetDistance.call(editor,point,frequency);
  assert.ok(Math.abs(distance-4)<1e-10,'packet width dragging inverts the displayed frequency span');
  if(layout===4)assert.ok(Math.abs(frequency-(8000+4*24.7*(1+.00437*8000)))<1e-9);
}

for (const depth of [0,.01,.1,.2,.3,.5,1])
  assert.ok(Math.abs(beatDepthValue(beatDepthPosition(depth))-depth)<1e-12);
assert.ok(beatDepthPosition(.3)>.5,'gentle beating gets most of the depth slider');

assert.equal(ringBeatRate(125,1.25,.5),1.25);
assert.equal(ringBeatRate(500,1.25,.5),2.5);
assert.equal(ringBeatRate(15000,80,1),80);

for (const rate of [0,.1,.5,1.25,3,10,80]) {
  assert.ok(Math.abs(beatRateValue(80,beatRatePosition(80,rate))-rate)<1e-10);
  assert.equal(hasPairedRing(3,rate),rate>0);
  assert.equal(hasPairedRing(4,rate),rate>0);
  assert.equal(hasPairedRing(2,rate),false);
}
assert.ok(beatRatePosition(80,3)>.5,'slow beating gets most of the slider travel');

for (let handles=0;handles<=32;handles++) {
  const points=Array.from({length:32},(_,i)=>({active:i<handles,allocation:1,width:.02}));
  for(const density of [0,.1,.45,.85,1]) {
    const result=packetAllocation(points,density,p=>p.width);
    assert.ok(result.count<=512 && result.count>=handles);
    assert.equal(result.count,handles+2*result.pairs.reduce((a,b)=>a+b,0));
    result.pairs.forEach((n,i)=>{if(i>=handles) assert.equal(n,0);});
    const paired=packetAllocation(points,density,p=>p.width,true);
    assert.ok(paired.count<=512);
    assert.equal(paired.count,2*handles+2*paired.pairs.reduce((a,b)=>a+b,0));
  }
}
const points=Array.from({length:12},()=>({active:true,allocation:1,width:.02}));
points[0].allocation=0; points[1].allocation=4;
const result=packetAllocation(points,1,p=>p.width);
assert.equal(result.count,512);
assert.equal(result.pairs[0],0);
assert.ok(result.pairs[1]>=3*result.pairs[2]);
console.log('Packet allocation: pool bounds, sparse/dense, inactive and weighted packets pass');
