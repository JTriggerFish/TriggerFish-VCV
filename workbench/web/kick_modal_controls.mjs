import {ModalEditor} from "./modal_editor.mjs";
import {mountModalTemplates} from "./modal_templates.mjs";

const Fields = ["frequency", "level", "centre", "edge"];
const key = (field, index) => `resonance_${field}_${index}`;

// Same editor as the metallic body; this resonator has no stochastic widths.
export class KickModalControls {
  constructor(owner) {
    this.owner = owner;
    const parent = document.getElementById("kick-modal-editor");
    parent.replaceChildren();
    this.editor = new ModalEditor(parent, {
      minimumFrequency:20, maximumFrequency:15000, frequencyScale:"log",
      minimumLevel:-72, maximumLevel:6, widthEnabled:false,
      points:()=>this.points(), globalTurbulence:()=>0, packetSpread:()=>0,
      replace:points=>this.replace(points),
      insert:(frequency, level)=>{
        const points=this.points(), index=points.findIndex(point=>!point.active);
        if(index<0)return null;
        points[index]={frequency,level,centre:1,edge:1,turbulence:0,active:true};
        this.replace(points); return index;
      },
      remove:index=>{
        const points=this.points(); points[index].level=-72; points[index].active=false;
        this.replace(points);
      },
      select:index=>this.inspect(index),
      readout:text=>{document.getElementById("kick-mode-readout").textContent=text;},
    });
    mountModalTemplates(document.getElementById("kick-modal-templates"), {
      capacity:16, minimumFrequency:20, maximumFrequency:15000,
      apply:generated=>this.replace(this.points().map((point,i)=>generated[i] ??
        {...point,level:-72,active:false})),
    });
    document.getElementById("kick-mode-edit").onclick=()=>this.editor.setTool("edit");
    document.getElementById("kick-mode-paint").onclick=()=>this.editor.setTool("paint");
    document.getElementById("kick-mode-clear").onclick=()=>this.replace(
      this.points().map(point=>({...point,level:-72,active:false})));
    this.inspect(null);
  }

  points() {
    return Array.from({length:16},(_,i)=>{
      const point=Object.fromEntries(Fields.map(field=>[
        field,this.owner.state.macros[this.owner.byKey.get(key(field,i)).index],
      ]));
      return {...point,active:point.level>-71.999,turbulence:0};
    });
  }

  replace(points) {
    points.forEach((point,i)=>Fields.forEach(field=>{
      const descriptor=this.owner.byKey.get(key(field,i));
      // Painting can create a point without spatial data. Start with unit
      // coupling rather than deriving new coefficients from its slot index.
      const value=point[field] ?? (field==="centre" || field==="edge" ? 1 : descriptor.defaultValue);
      this.owner.state.macros[descriptor.index]=Math.max(descriptor.minimum,
        Math.min(descriptor.maximum, value));
    }));
    this.editor.refresh();
    this.inspect(this.editor.selected);
    this.owner.onChange("resonance_modes");
  }

  inspect(index) {
    const parent=document.getElementById("kick-mode-selection");
    parent.replaceChildren();
    const title=document.createElement("b");
    title.textContent=index===null ? "No mode selected" : `Mode ${index+1}`;
    parent.append(title);
    const names=["Frequency","Relative prominence","Centre strike coupling","Edge strike coupling"];
    Fields.forEach((field,j)=>{
      this.owner.addSlider("kick-mode-selection",key(field,index??0));
      const row=parent.lastElementChild;
      row.querySelector("span").textContent=names[j];
      const input=row.querySelector("input");
      input.disabled=index===null;
      const handler=input.oninput;
      input.oninput=()=>{handler();this.editor.refresh();};
      const reset=input.ondblclick;
      input.ondblclick=event=>{reset(event);this.editor.refresh();};
      if(field==="centre" || field==="edge")
        row.dataset.tooltip="Signed strike coupling at this location endpoint. Negative values invert polarity; the body normalizes the combined excitation energy.";
      if(index===null)row.querySelector("output").textContent="—";
    });
  }

  destroy() { this.editor.resizeObserver.disconnect(); }
}
