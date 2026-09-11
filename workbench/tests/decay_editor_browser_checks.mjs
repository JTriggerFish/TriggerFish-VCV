// Real pointer checks in the caller's disposable, silent workbench tab.
import assert from "node:assert/strict";

export async function checkDecayInteractions(call, evaluate) {
  await evaluate(`(async()=>{
    document.getElementById('decay-editor').scrollIntoView({block:'center'});
    await new Promise(requestAnimationFrame);
  })()`);
  const locate = selector => evaluate(`(()=>{
    const box=document.querySelector(${JSON.stringify(selector)}).getBoundingClientRect();
    return {x:box.x+box.width/2,y:box.y+box.height/2};
  })()`);
  const seconds = () => evaluate(`parseFloat(document.querySelector(
    '#decay-selection [data-fit-key^=body_decay_seconds_] output').textContent)`);
  const drag = async (selector, dy, modifiers=0) => {
    const p=await locate(selector);
    await call('Input.dispatchMouseEvent',{type:'mousePressed',...p,button:'left',clickCount:1,modifiers});
    await call('Input.dispatchMouseEvent',{type:'mouseMoved',x:p.x,y:p.y+dy,buttons:1,modifiers});
    await call('Input.dispatchMouseEvent',{type:'mouseReleased',x:p.x,y:p.y+dy,button:'left',clickCount:1,modifiers});
    return p;
  };
  const selector='#decay-editor [data-slot="0"]';
  const initial=await seconds();
  const start=await drag(selector,-10);
  const raised=await seconds(), moved=await locate(selector);
  assert.ok(raised>initial,'Dragging up must lengthen T60');
  assert.ok(Math.abs(moved.y-start.y+10)<.2,'Knot must follow pointer on softened-log scale');
  await drag(selector,10);
  assert.ok(Math.abs(await seconds()-initial)<.02,'Reverse drag must restore T60');
  const fineStart=await locate(selector);
  await drag(selector,-10,8); // CDP Shift modifier.
  const fineEnd=await locate(selector);
  assert.ok(Math.abs(fineEnd.y-fineStart.y+1)<.2,'Shift must give tenfold finer travel');
  await drag(selector,10,8);

  const low=await seconds();
  const highHandle='#decay-editor [data-slot="7"]';
  await drag(highHandle,0);
  const high=await seconds();
  await drag('#decay-editor .decay-all-handle',-5);
  const highAfter=await seconds();
  await drag(selector,0);
  const lowAfter=await seconds();
  assert.ok(Math.abs(lowAfter/low-highAfter/high)<.02,'Diamond must preserve decay ratios');
  assert.ok(lowAfter>low && highAfter>high);

  await evaluate(`(async()=>{
    const amount=document.querySelector('[data-fit-key=field_motion_depth] input');
    amount.value=.1;amount.dispatchEvent(new Event('input'));
    const value=parseFloat(amount.nextElementSibling.textContent);
    if(Math.abs(value-.03)>1e-6)throw Error('Shimmer did not use squared mapping');
    await document.getElementById('instrument-calibration').onchange();
    document.querySelector('aside').scrollTop=0;
  })()`);
  console.log('T60 pointer tracking, reverse/fine drags, ratio-preserving shift and shimmer scaling pass');
}
