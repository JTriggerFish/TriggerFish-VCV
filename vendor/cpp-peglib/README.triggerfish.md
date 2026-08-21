# cpp-peglib vendoring record

- Upstream: https://github.com/yhirose/cpp-peglib
- Release: `v1.9.1`
- Commit: `527422aa38ce78c4214b21843e97584bc2db0150`
- `peglib.h` SHA-256: `076A67042C0F703B9B58F2077C4A68AAF68E46FE2D6DA4D8A134616E2A6EE433`
- License: MIT; the unmodified upstream `LICENSE` is stored beside the header.

The dependency is vendored because VCV Rack plugins must not depend on a
parser being installed on the user's system. Update it only by choosing a
specific upstream release, replacing both upstream files, recording the new
commit and hashes here, and running the sequencer parser tests.
