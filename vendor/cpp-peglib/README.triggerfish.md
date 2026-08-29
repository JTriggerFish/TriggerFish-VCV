# cpp-peglib vendoring record

- Upstream: https://github.com/yhirose/cpp-peglib
- Release: `v1.9.1`
- Commit: `527422aa38ce78c4214b21843e97584bc2db0150`
- `peglib.h` upstream SHA-256: `076A67042C0F703B9B58F2077C4A68AAF68E46FE2D6DA4D8A134616E2A6EE433`
- `peglib.h` patched SHA-256: `D621E2FE3CAE59F8DEE64A4E9A9FCCDF09831E80808E9B427CE448EAE388CB45`
- License: MIT; the unmodified upstream `LICENSE` is stored beside the header.

The vendored header includes cpp-peglib's former header-only `any`
implementation, backported from upstream release `v0.1.14` (commit
`b92da07beddd286cc16ef3620793e297a5f17a6c`). It is selected for macOS
deployment targets older than 10.13 because `std::any` in the macOS 12.3 SDK is
unavailable at VCV Rack's macOS 10.9 deployment target. CI exercises both the
standard-library and compatibility implementations.

The dependency is vendored because VCV Rack plugins must not depend on a
parser being installed on the user's system. Update it only by choosing a
specific upstream release, replacing both upstream files, reapplying and
documenting the compatibility patch, recording the new commit and hashes here,
and running the sequencer parser tests with both `any` implementations.
