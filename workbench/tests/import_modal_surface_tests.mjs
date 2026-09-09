import assert from 'node:assert/strict';
import {importModalSurface,importOutputEq} from '../web/import_modal_surface.mjs';
import {turbulenceIntensity} from '../web/turbulence_profile.mjs';
const descriptors=[{key:'field_wander_hz'}];
const fit=p=>({instrument:{recipe:'metal.cymbal.v1',nodes:[{id:'body',parameters:p}]}});
for (const centre of [1,880,4000,15000]) for (const slope of [-1,0,.4,1]) {
  const original=fit({field_turbulence:4,field_turbulence_centre:centre,
    field_turbulence_slope:slope,field_drift_depth:0,field_drift_rate:8});
  const converted=importModalSurface(original,descriptors).instrument.nodes[0].parameters;
  assert.equal(original.instrument.nodes[0].parameters.field_turbulence_centre,centre);
  assert.ok(!Object.hasOwn(converted,'field_turbulence_centre'));
  assert.equal(converted.field_wander_hz,0);
  for(const f of [20,125,1000,15000]) {
    const old=4*(f/centre)**slope;
    const now=turbulenceIntensity(f,converted.field_turbulence,slope,1000,1,true);
    assert.ok(Math.abs(now-old)<1e-9*Math.max(old,1));
  }
}
assert.throws(()=>importModalSurface(fit({field_drift_depth:.1}),descriptors),/cannot be converted/);
assert.throws(()=>importModalSurface(fit({field_turbulence:NaN,field_turbulence_centre:500,
  field_turbulence_slope:0}),descriptors),/Invalid legacy/);
console.log('modal surface imports preserve curves, reject lossy drift conversion');

const eqDescriptors=[{key:'output_eq_enabled'}];
const eqFit=p=>({instrument:{recipe:'metal.cymbal.v1',nodes:[{id:'observation',parameters:p}]}});
const oldEq={direct_gain:.18,field_gain:1.7,
  direct_radiation_enabled:0,direct_low_cut:80,direct_colour_frequency:1200,
  direct_colour_gain:-2,direct_high_cut:5000,
  body_radiation_enabled:1,body_low_cut:25,body_colour_frequency:7200,
  body_colour_gain:2,body_high_cut:19000};
const original=eqFit(oldEq);
const converted=importOutputEq(original,eqDescriptors);
assert.deepEqual(converted.instrument.nodes[0].parameters,{
  direct_gain:.18,field_gain:1.7,output_eq_enabled:1,output_low_cut:25,
  output_colour_frequency:7200,output_colour_gain:2,output_high_cut:19000});
assert.deepEqual(original,eqFit(oldEq));
assert.equal(importOutputEq(converted,eqDescriptors),converted);
assert.throws(()=>importOutputEq(eqFit({...oldEq,body_low_cut:NaN}),eqDescriptors),/invalid old EQ/);
assert.throws(()=>importOutputEq(eqFit({...oldEq,output_eq_enabled:1}),eqDescriptors),/mixes per-path/);
assert.throws(()=>importOutputEq(eqFit({body_low_cut:25}),eqDescriptors),/Missing/);
console.log('shared EQ import removes old controls, preserves levels and body curve');
