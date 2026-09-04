# Reference listening targets

The development workbench exposes one deliberately simple reference pairing
per instrument. Each entry selects a compiled recipe, one documented starting
state, one standard articulation at a middle velocity, and the exact local
reference cell. Obsolete fitted metallic presets are deliberately not migrated
or retained.

| Toolbar entry | Recipe | Reference target | Parameter state |
| --- | --- | --- | --- |
| Crash — medium edge | `metal.cymbal.v1` | edge, velocity 72, take 1 | current factory defaults |
| Snare — medium standard hit | `drum.snare.v1` | main, velocity 82, take 1 | current factory defaults |
| Acoustic kick — medium centre | `drum.membrane.v1` | centre, velocity 64, take 1 | -12 dB; explicit beater, FM/body, membrane, and observation routes |
| Gong — representative mallet | `metal.cymbal.v1` | mallet, velocity 96, take 3 | measured `gong-v1` start |
| Ride — medium bow | `metal.cymbal.v1` | bow, velocity 82, take 1 | current factory defaults |
| Hi-hat — medium half-open | `metal.cymbal.v1` | half-open, velocity 96, take 1 | current factory defaults plus passive pedal constraint |

The highlighted **Reference targets** chooser is kept in the persistent top toolbar,
next to the recipe selector. Loading an entry performs one atomic user action:

1. activate the required compiled recipe;
2. restore its documented current starting state;
3. select and load the named reference cell;
4. copy the cell's strength, location, implement, hardness, spread, passive
   constraint/pedal position, and seed to
   the performance controls; and
5. preserve the collection-level reference monitoring calibration.

None of these entries is an accepted synthesis calibration. In particular, a
parameter vector produced by an optimizer is not promoted merely because it
improves a relative or aggregate score.

Each factory model level can be checked against its named reference without changing
either waveform. It is not a per-hit normalizer and it is not evidence that the
other velocity layers match. The diagnostic command
`./dev.ps1 audit-workbench-starts` reports raw whole-response energy, peak, and
causal-region band errors for the exact served patch. Those numbers can reject
an obviously broken start; only listening can accept one.

New fitted states must be saved explicitly against the current one-state modal
topology. Development snapshots from retired schemas are rejected rather than
silently translated.

Velocity is never passed through saturation, a soft knee, or a fitted amplitude
exponent. Every instrument's base event amplitude is linear. Individual bright
sources such as kick click or cymbal micro-contact may grow faster as an
explicit timbral response, while velocity-dependent pitch and spectral
excitation remain available because they change articulation rather than
flattening level.

The ordinary reference selectors remain available for neighbouring cells.
Corpus monitor gains are fixed once per source collection and affect reference
and synthesis playback equally. They compensate large collection-level
monitoring differences; they are not fitted synthesis parameters and never
normalize individual files.
The model-level control is attenuation-only (`-60..0 dB`). If matching would
need positive gain, the workbench stops at 0 dB and reports that the starting
point itself is under-levelled instead of hiding the defect.
Double-clicking the model-level control remains an explicit one-cell comparison
aid. Loading a reference target or moving between velocity cells does not invoke
that operation, because per-cell normalization would erase the velocity curve.
The current private crash view is intentionally limited to five velocities,
five articulations, and one repeat. Snare, kick, gong, ride, and hi-hat use similarly
small allow-listed views so the browser stays useful rather than becoming a
sample-library explorer.

The hi-hat is deliberately a parameterization of `metal.cymbal.v1`, not a
parallel synthesis engine. Its pedal is the shared plate constraint: open and
partly closed states change modal loss continuously, while the same contact and
stochastic modal field remain in use. A constraint already present before a hit
is installed as an initial condition; changes during a sound are smoothed and
strictly dissipative. A future foot-chick source will be a contact gesture into
this same model, not a delayed playback event.

Sample data and machine paths are never serialized into the repository. The
catalog stores only neutral corpus aliases, public event mappings, and relative
HTTP identifiers resolved by the local development server. The browser
workbench, WebAssembly engine, and these calibration records remain optional
development targets and are not dependencies of Rack or release builds.
