# Ride cymbal reference corpus

This note records a deliberately multi-source corpus.  No one library is the
fitting truth: model trends must survive different cymbals, recording chains,
articulations, and sample-library mappings.

## Primary fitting references

| Source | Coverage | Licence / intended use | Role |
| --- | --- | --- | --- |
| [Versilian Community Sample Library (VCSL)](https://github.com/sgossner/VCSL) | Two suspended cymbals; stick and mallet dynamics; bell hits | CC0; its authors explicitly permit commercial software | Exciter hardness, bell/body contrast, and velocity response |
| [Salamander Drumkit](https://github.com/endolith/Salamander-Drumkit) | Two 20-inch rides; 51 bow, bell, edge/crash, choke, pp/mp/ff and repeated-take WAVs | CC BY-SA | Ride decay, variation, velocity, and articulation baseline |
| [29kSamplesDrumsDataset](https://zenodo.org/records/4958592) | Ride hits from two cymbal sets, recorded with eight microphones | CC BY 4.0 | Population statistics and microphone-robust descriptors |
| [University of Iowa MIS](https://theremin.music.uiowa.edu/MIScymbals.html) | One 21-inch ride; bell/normal/shoulder at pp/mf/ff, mallet and roll; anechoic, documented chain | Described as freely available, but no explicit content licence was found | Selected v1 object: bow/bell MF/FF structural, damping, location, and velocity fit |
| [Bitwig Acoustic Drums and Percussion](https://www.bitwig.com/sound-content/acoustic-drums-and-percussion-127/) | Four contrasting 20-inch rides, plus 18-inch crash/ride and sizzle ride; seven variation groups with four dynamics each | Installed licensed content; source audio remains local and is not redistributed | Held-out multi-model velocity and population validation; no inferred bell labels |

The local, ignored working copies live below
`build/cymbal-calibration/references`. Do not commit third-party audio.

## Frozen v1 primary object

`data/ride-calibration/iowa-21ride-mf-ff-v1.json` selects exactly four source
AIFFs from one Iowa 21-inch ride: normal/bow MF and FF, and bell MF and FF. The
source hashes, source labels, model-strength mapping, fixed channel transform,
and exclusions are tracked. The replacement dataset loader must verify those
hashes and produce 48 kHz mono float WAVs using
`mid = 0.5 * (left + right)`, with no level normalization. Derived hashes and
qualification measurements belong in an ignored generated manifest at
`build/cymbal-calibration/fit-targets/iowa-21ride-mf-ff-v1/manifest.json`.

The rebuilt set passes the selection gates:

- bow and bell MF-to-FF early levels rise by 22.50 and 19.79 dB;
- the two location-specific velocity gains differ by only 2.71 dB;
- bell sustained high-mid share exceeds bow by 7.38 dB at MF and 5.27 dB at FF;
- every cell has at least 10.06 seconds after detected onset; and
- the original L/R correlation is 0.76--0.82, channel imbalance is at most
  0.24 dB, and the declared mid downmix loses at most 0.56 dB, so mono
  conversion is not cancelling a strong out-of-phase component; and
- all sources have one known physical object, recording chain, implement, and
  labelled strike region.

The earlier derived Iowa WAVs are not fit inputs: FF had accidentally been
converted to mono 16-bit while MF remained stereo float. The v1 builder starts
from the consistently stereo 32-bit/44.1 kHz AIFF sources and makes the mono
policy explicit.

## Perceptual cross-checks

- [Paiste](https://www.paiste.com/en/products/models/108-ride) publishes bell,
  body, edge, and pattern MP3s for a large range of rides.
  It includes unusually useful controlled size families: Dark Energy Mark I and
  II at 20/21/22 inches, Formula 602 rides at 20/22/24 inches, 2002 rides at
  20/22/24 inches, and several flat/dry/power contrasts. Paiste labels the files
  as intended for authorized use only.
- [Zildjian](https://zildjian.com/products/20-k-custom-ride) product pages
  commonly provide bell, crescendo, and ride-pattern
  clips. Useful contrasts include the 20-inch K Custom Ride and 21-inch A Sweet
  Ride. These are product demonstrations rather than controlled measurements.
- [Sabian](https://sabian.com/product/22012-20-inch-aa-medium-ride/) product
  pages commonly provide separate bow, bell, and performance
  clips and publish size, weight, alloy, and timbre metadata. The 21-inch HH Raw
  Bell Dry Ride and 22-inch HH Power Bell Ride form a useful dry/power contrast.
- Meinl and Istanbul publish useful product demonstrations and construction
  metadata, but the current pages are less convenient for isolated, labelled
  one-shots.

Manufacturer MP3s may be normalized, compressed, processed, and recorded through
different signal chains. Use them to rank spectral envelope, articulation
contrast, and broad decay tendencies. Do not fit absolute level or fine modal
frequencies across brands from these clips.

## Best paid complement

[Sonomar Collection: Cymbals](https://sonomar.ca/en/sonomar-collection-cymbals/)
is the most useful conventional purchase found so far.  It contains 1,386
24-bit/192-kHz files with embedded metadata, made from 450 unique recordings
through three stereo microphone pairs.  The team varied bell/rim/bow location,
wood-stick and mallet implement, speed, angle, and intensity.  Its rides include
a 20-inch Zildjian K Custom Medium, 21-inch Paiste Dark Energy Mark I, and
24-inch Paiste Giant Beat.  That is unusually good coverage for location,
implement, size, and microphone robustness.  The current local/open/Bitwig
corpus is already broad enough to build the analysis pipeline, so purchase this
for a later independent fit/validation expansion rather than making it the
first or only target.

## Installed proprietary libraries

Superior Drummer 3 Core contains these stick rides:

- 20-inch Zildjian Kerope
- 22-inch Zildjian K Constantinople Medium Thin Low
- 20-inch Istanbul Agop Turk
- 24-inch Paiste 2002
- 22-inch Istanbul Agop 30th Anniversary

The Kerope also has felt-mallet, brush, and rod recordings. Its ride mapping
exposes bow tip (MIDI 51), bow shank (29), bell tip (30), bell shank (53), and
edge (59). EZdrummer 3 adds 21/22/24-inch Paiste Dark Energy, Formula 602 and
2002 variants, a 24-inch 2002 Swish, a 24-inch Zildjian A Trans Stamp, and a
22-inch Istanbul Agop 30th Anniversary variant. Its common mapping is bow 51,
bell 53, and edge 59.

`tools/build_toontrack_ride_sweep.py` generates deterministic MIDI and an onset
JSON manifest for private listening tests. Toontrack's current EULA restricts
using its product sounds to develop distributed software and restricts uses
outside musical works without written permission. Do not use these renders as a
reproducible fitting dataset unless Toontrack grants written permission.

Bitwig's installed Acoustic Drums and Percussion package contains:

- 18-inch K Custom crash/ride
- 20-inch A Custom ride
- 20-inch K ride
- 20-inch K Custom ride
- 20-inch K Custom Dry ride
- K Custom sizzle ride

Each conventional 20-inch ride has seven groups (A-G) with four clear dynamic
steps per group. Measurements across all four rides show that the letter-group
contrasts are not stable radial-location contrasts: no group has a consistent
bell-versus-bow signature. Treat A-G as take/strike variation populations, not
as bell/bow/edge labels. Their four monotonic layers remain useful for held-out
velocity and density validation. Do not copy the source recordings into the
repository or distribute them. The fit must also hold across independent
CC0/CC-BY/CC-BY-SA references so that it does not reproduce this library's
recording chain.

## Capture recommendation

Before buying another conventional drum library, ask the vendor this exact
question: may isolated samples be used privately to measure and tune the
parameters of a commercial physical/perceptual cymbal model, provided no source
audio or recoverable derivative samples are distributed? A normal music-use
licence is not enough.

If no vendor grants this explicitly, commission a compact recording session. A
high-value matrix is four contrasting rides (roughly 18/20/22/24 inches), six
locations/implements (bell tip and shank, inner/middle/outer bow tip, edge
shank), eight controlled velocities, and four repeats. Record a close channel
and a neutral overhead at 24-bit/96 kHz, with processing off and enough time for
the full tail. This is smaller than a premium drum library, has unambiguous
metadata, and gives us permanent model-development rights.
