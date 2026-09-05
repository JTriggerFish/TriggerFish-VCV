# Modular percussion patches

## Purpose

A TriggerFish percussion instrument is a versioned JSON patch that connects
reusable contact, body, interaction, observation, and output modules. The JSON
is the durable design format shared by the browser workbench, fitting tools,
factory instruments, and Rack serialization. It is not the audio-thread
execution representation.

The first three production recipes are the metallic plate, compact FM kick,
and shared-state membrane drum. All expose named patch parameters through one
recipe-aware C/WebAssembly API;
the earlier crash-only API remains only as a compatibility wrapper.

## Signal-flow boundary

The visible instrument graph is a directed acyclic graph:

```text
events -> exciters -> transformations -> bodies -> interactions
                \-----------------------------> observation -> output
```

Feedback is owned by a component with a tested recurrence. A dispersion loop,
feedback-delay body, or locally coupled modal field is therefore one node, not
a cycle in the patch graph. This keeps scheduling deterministic and prevents a
patch editor from creating an unqualified zero-delay feedback loop.

Modules have typed, named input and output ports. The membrane recipe adds the
first scalar-control port between its strike-history/tension envelope and modal
body. That envelope is not presented as measured body energy; event,
body-motion, and normalized-energy ports are added only with concrete
consumers. Port types, component versions,
parameter descriptors, and presentation metadata belong to one module
registry rather than being repeated in every patch.

## JSON contract

The initial schema identifier is `triggerfish.percussion.patch/v1`. A patch
contains:

- a stable patch ID and display name;
- a minimum compatible engine version;
- a registered production recipe;
- versioned module instances with stable IDs and named parameters;
- typed connections and named outputs;
- reduced performance controls and their mappings; and
- optional editor layout metadata ignored by DSP.

Parameter names include units where ambiguity is possible. Positional arrays
are transport details and must not become the durable patch format. Unknown
module versions, parameter names, ports, cycles, duplicate IDs, non-finite
values, and out-of-range values are rejected before a voice is prepared.
Every connection contains one explicit Boolean `enabled` state and no numeric
gain. Audible levels are named parameters owned by modules, so the UI, JSON,
and C++ parameter surfaces cannot conceal an additional edge coefficient.
Deterministic expansion of an explicitly labelled macro is different: `Size
meta` rewrites visible parameters, while event controls such as `Implement`
derive per-hit contact projections without creating another saved control.

Fitting snapshots embed this patch alongside their reference identity, event,
and analysis settings. The patch is therefore reproducible without relying on
the current value of a factory preset or an external file.

## Workbench execution

The browser first validates the patch, then matches it to one of a bounded set
of registered recipes. A recipe owns a statically ordered schedule of typed C++
calls; it is not a general graph interpreter. The metallic-plate recipe
compiles three optional audio connections to prepared switches at those call sites.
The unified kick recipe similarly compiles three observation switches around
one contact, one swept-sine thump and a persistent membrane resonator. See
[kick architecture](TfPercussion-kick-architecture.md). The membrane recipe compiles
five optional switches around contact and correlated-FM exciters. Four visible
source-to-bus level parameters feed its two fixed mixers, followed by a
normalized and energy-bounded persistent 16-mode membrane, a
strike-history/tension envelope, a
two-source observation, and a selectable output EQ. This makes every routing
state explicit while keeping the sample loop free of JSON parsing, allocation,
virtual dispatch, and topology discovery.

The workbench may enable or disable each recipe's optional source connections
live. It
rejects a routing state with no complete path to the output, and required
mix/observation/output connections are locked. Module replacement and arbitrary
new connections remain unavailable until another registered recipe declares
and implements them. A structural edit prepares a replacement voice away from
the audio callback and
swaps it at a safe boundary. Expensive membrane modal coefficients are carried
in the prepared recipe blob, and fixed-capacity observation delays keep
prepared installation allocation-free. Parameters declare whether they are structural,
sampled at the next hit, or continuously smoothed.

The compact routing overview remains visible above the complete parameter
panel. Double-clicking it opens the larger design view. Module role colours are
shared by routing nodes and parameter-section markers, but text labels remain
the primary identifier. All modules present in the patch keep their controls
visible; graph selection only locates or highlights a section.

## Rack execution

Factory Rack instruments use registered, statically compiled topology recipes
with the JSON patch embedded as prepared parameter data. The
sample loop remains direct function calls: it performs no JSON parsing,
topology discovery, allocation, or general graph traversal.

A Rack module exposes only the patch's `performanceControls`. One performance
control may map through bounded curves to several module parameters—for
example, decay can scale body T60 and wire persistence by different amounts.
Persistent knobs may be overridden by CV, while hit strength remains an event
input. Rack project serialization embeds the resolved instrument RTTI-free
patch data so later factory-preset changes do not alter an existing project.

The development-only WebAssembly, JavaScript UI, reference corpus, and analysis
tools remain separate build targets and are not dependencies of normal Rack or
official plugin builds.

## Recipe family

The intended initial recipes are:

| Recipe | Construction |
| --- | --- |
| Metallic plate | contact, stochastic modal field with intrinsic energy cascade, observation |
| Membrane | contact, modal or feedback body, optional tension state, observation |
| Snare | membrane plus body-driven wire interaction |
| Compact kick | two overlap-safe correlated-FM burst branches |
| Acoustic kick | reduced FM/contact source feeding an optional head/cavity body |
| Hi-hat | metallic body with pedal-controlled passive loss and driven plate contacts |

These are convenient starting structures, not hard instrument categories. The
same components may form gongs, toms, rattles, prepared drums, or synthetic
hybrids as long as their port contracts and energy invariants remain valid.

## Implementation sequence

1. **Complete:** qualify the patch schema, validator, migration, compact routing
   overview, and executable metallic-plate routing variants.
2. **Complete:** move shared module descriptors and recipe matching behind the
   Wasm boundary, retain statically ordered schedules, and prove the boundary
   with a second complete compact-kick recipe.
3. **Complete:** add reusable fixed mixers, transient voice pooling, and
   selectable bypass/radiation/multiband observation EQ.
4. **Complete:** implement the shared-state modal membrane recipe and expose
   tom and acoustic-kick starting patches in the workbench.
5. **Partial:** the strike-history tension envelope is implemented; add the body-driven
   snare-wire interaction.
6. Generate reduced Rack controls from performance mappings and embed factory
   patches in the instrument modules.
