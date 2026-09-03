# Reference calibration starting points

The development workbench exposes one deliberately simple reference pairing
per instrument. Each entry selects a compiled recipe, a complete visible
parameter state, one standard articulation at a middle velocity, and the exact
local reference cell. These are starting points for listening work, not claims
that every instrument or dynamic layer is fully fitted.

| Toolbar entry | Recipe | Reference target | Parameter state |
| --- | --- | --- | --- |
| Crash — medium edge | `metal.cymbal.v1` | edge, velocity 72, take 1 | metallic factory starting point |
| Snare — medium centre | `drum.snare.v1` | main, velocity 82, take 1 | initial acoustic-snare fit |
| Kick — medium centre | `drum.kick-fm.v1` | centre, velocity 64, take 1 | compact FM kick starting point |
| Gong — representative mallet | `metal.cymbal.v1` | mallet, velocity 96, take 3 | large-metal size endpoint |
| Ride — medium bow | `metal.cymbal.v1` | bow, velocity 82, take 1 | ride-specific modal/decay starting point |

The highlighted **Calibrations** chooser is kept in the persistent top toolbar,
next to the recipe selector. Loading an entry performs one atomic user action:

1. activate the required compiled recipe;
2. restore its parameter starting point;
3. select and load the named reference cell;
4. copy the cell's strength, location, implement, hardness, spread, and seed to
   the performance controls; and
5. request the standard model-to-reference level match.

The ordinary reference selectors remain available for neighbouring cells.
Corpus trims are fixed per source and apply equally to reference and synthesis
in the protected listening path; they never peak-normalize individual files.
The current private crash view is intentionally limited to five velocities,
five articulations, and one repeat. Snare, kick, gong, and ride use similarly
small allow-listed views so the browser stays useful rather than becoming a
sample-library explorer.

Sample data and machine paths are never serialized into the repository. The
catalog stores only neutral corpus aliases, public event mappings, and relative
HTTP identifiers resolved by the local development server. The browser
workbench, WebAssembly engine, and these calibration records remain optional
development targets and are not dependencies of Rack or release builds.
