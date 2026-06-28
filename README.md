# Quotient

This is a work in progress!

A dynamic EQ built on numerically stable, trapezoidally-integrated state-variable filters and selectable models of classic gain-reduction behavior.

## Overview

Quotient pairs per-band dynamic processing with a choice of phase responses:

- **State-variable filters** discretized via the topology-preserving transform (trapezoidal integration). Unconditionally stable, with low coefficient sensitivity across the full frequency range.
- **Selectable gain-reduction models** reproducing the response of well-known dynamics processors, applied per band.
- **Minimum-phase IIR** for low-latency processing with no pre-ringing.
- **Linear-phase FIR** via partitioned (block) convolution using FFT/IFFT, for phase-coherent processing where latency is acceptable.

## Status

Quotient is primarily a vehicle for rigorously exercising applied mathematics.

## Build

Quotient builds with CMake and JUCE. All third-party code lives in `libs/` as git submodules.

### Quick start
```sh
git clone --recurse-submodules <repo-url> Quotient
cd Quotient
./scripts/bash/setup.sh      # installs deps + initialises submodules (one-time)
./scripts/bash/build.sh      # configures and compiles (Debug by default)
```
Run `./scripts/bash/setup.sh` if you've already cloned without submodules.

### Prerequisites
`setup.sh` installs these for you on macOS via Homebrew:
- A C++23 toolchain (Xcode Command Line Tools on macOS)
- CMake ≥ 3.23 and Ninja
- Boost (`find_package(Boost REQUIRED)`)
- Python 3 (`find_package(Python3 REQUIRED)`)
- Doxygen (optional, for the `docs` target)

### build.sh options
- `./scripts/bash/build.sh --release` — optimised build
- `./scripts/bash/build.sh --target Quotient_VST3` — build a single format/target
- `./scripts/bash/build.sh --clean` — wipe the build directory and reconfigure
- `./scripts/bash/build.sh --help` — full option list

Built plugins are written to `build/Quotient_artefacts/<Config>/` as **Standalone**, **VST3**, and (on macOS) **AU**.

## Feedback

Notes and corrections are welcome! Please reach out.
