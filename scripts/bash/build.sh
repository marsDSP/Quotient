#!/usr/bin/env bash
#
# build.sh — configure + compile the Quotient JUCE plugin.
#
# Designed so that a fresh clone builds with a single command, and so that the
# most common CMake/JUCE pitfalls are either avoided or explained:
#   * missing/uninitialised git submodules (libs/JUCE …)
#   * no Ninja on PATH (falls back to Unix Makefiles)
#   * stale CMake cache / generator mismatch (auto-reconfigures)
#   * CMake 4 dropping pre-3.5 policy support in vendored libs
#   * Boost not found (adds the Homebrew prefix to CMAKE_PREFIX_PATH)
#
# Usage:
#   ./scripts/bash/build.sh [options]
#
# Options:
#   --debug              Configure a Debug build (default).
#   --release            Configure a Release build.
#   -t, --target <name>  Build a specific target (e.g. Quotient_VST3,
#                        Quotient_AU, Quotient_Standalone, Quotient_All,
#                        AudioPluginHost, RunAudioPluginHost, docs).
#   -b, --build-dir <d>  Build directory (default: build).
#   -G, --generator <g>  CMake generator override (default: auto-detect).
#   -j, --jobs <n>       Parallel build jobs (default: CPU count).
#   -c, --clean          Remove the build directory before configuring.
#   -h, --help           Show this help.
#
set -euo pipefail

# --------------------------------------------------------------------------
# Paths / defaults
# --------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# This script lives in <repo>/scripts/bash, so the project root (where
# CMakeLists.txt, libs/ and .gitmodules live) is two directories up.
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$ROOT"

BUILD_DIR="build"
BUILD_TYPE="Debug"
TARGET=""
GENERATOR=""
JOBS=""
CLEAN=0

# Public upstream JUCE fallback
JUCE_FALLBACK_URL="https://github.com/juce-framework/JUCE.git"

# --------------------------------------------------------------------------
# Logging helpers
# --------------------------------------------------------------------------
if [[ -t 1 ]]; then
  C_BLUE="\033[1;34m"; C_GREEN="\033[1;32m"; C_YELLOW="\033[1;33m"
  C_RED="\033[1;31m"; C_RESET="\033[0m"
else
  C_BLUE=""; C_GREEN=""; C_YELLOW=""; C_RED=""; C_RESET=""
fi
info() { printf "${C_BLUE}==>${C_RESET} %s\n" "$*"; }
ok()   { printf "${C_GREEN}  ✓${C_RESET} %s\n" "$*"; }
warn() { printf "${C_YELLOW}  ! %s${C_RESET}\n" "$*" >&2; }
die()  { printf "${C_RED}  ✗ %s${C_RESET}\n" "$*" >&2; exit 1; }

usage() {
  cat <<'EOF'
build.sh — configure + compile the Quotient JUCE plugin.

Builds a fresh clone with a single command and handles the most common
CMake/JUCE pitfalls automatically:
  * uninitialised git submodules (auto-runs submodule update / JUCE clone)
  * no Ninja on PATH (falls back to Unix Makefiles)
  * stale CMake cache / generator mismatch (auto-reconfigures)
  * CMake 4 dropping pre-3.5 policy support in vendored libs
  * Boost not found (adds the Homebrew prefix to CMAKE_PREFIX_PATH)

Usage:
  ./scripts/bash/build.sh [options]

Options:
  --debug              Configure a Debug build (default).
  --release            Configure a Release build.
  -t, --target <name>  Build a specific target (e.g. Quotient_VST3,
                       Quotient_AU, Quotient_Standalone, Quotient_All,
                       AudioPluginHost, RunAudioPluginHost, docs).
  -b, --build-dir <d>  Build directory (default: build).
  -G, --generator <g>  CMake generator override (default: auto-detect).
  -j, --jobs <n>       Parallel build jobs (default: CPU count).
  -c, --clean          Remove the build directory before configuring.
  -h, --help           Show this help.
EOF
  exit 0
}

# --------------------------------------------------------------------------
# Argument parsing
# --------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --debug)         BUILD_TYPE="Debug" ;;
    --release)       BUILD_TYPE="Release" ;;
    -t|--target)     TARGET="${2:?--target requires a value}"; shift ;;
    -b|--build-dir)  BUILD_DIR="${2:?--build-dir requires a value}"; shift ;;
    -G|--generator)  GENERATOR="${2:?--generator requires a value}"; shift ;;
    -j|--jobs)       JOBS="${2:?--jobs requires a value}"; shift ;;
    -c|--clean)      CLEAN=1 ;;
    -h|--help)       usage ;;
    *) die "Unknown option: $1 (try --help)" ;;
  esac
  shift
done

# --------------------------------------------------------------------------
# Preflight: cmake + submodules
# --------------------------------------------------------------------------
command -v cmake >/dev/null 2>&1 || die "cmake not found. Run ./setup.sh first."
info "Using $(cmake --version | head -1)"

# JUCE is required by add_subdirectory(libs/JUCE). Auto-recover a checkout that
# was cloned without --recurse-submodules.
if [[ ! -f "$ROOT/libs/JUCE/CMakeLists.txt" ]]; then
  warn "libs/JUCE is missing — initialising git submodules…"
  if [[ -f "$ROOT/.gitmodules" ]]; then
    git submodule sync --recursive || true
    git submodule update --init --recursive || warn "submodule update reported errors."
  fi
  if [[ ! -f "$ROOT/libs/JUCE/CMakeLists.txt" ]]; then
    warn "Cloning JUCE directly from $JUCE_FALLBACK_URL"
    rm -rf "$ROOT/libs/JUCE"
    git clone --depth 1 "$JUCE_FALLBACK_URL" "$ROOT/libs/JUCE" \
      || die "Could not obtain JUCE. Run ./setup.sh and check your network/credentials."
  fi
fi

# --------------------------------------------------------------------------
# Generator selection (prefer Ninja, fall back to Unix Makefiles)
# --------------------------------------------------------------------------
if [[ -z "$GENERATOR" ]]; then
  if command -v ninja >/dev/null 2>&1; then
    GENERATOR="Ninja"
  else
    GENERATOR="Unix Makefiles"
    warn "ninja not found — using 'Unix Makefiles' (install ninja for faster builds: brew install ninja)."
  fi
fi

# Parallel jobs
if [[ -z "$JOBS" ]]; then
  JOBS="$( (command -v sysctl >/dev/null 2>&1 && sysctl -n hw.ncpu) \
          || (command -v nproc >/dev/null 2>&1 && nproc) \
          || echo 4 )"
fi

# Help find Homebrew-installed Boost (find_package(Boost REQUIRED)).
PREFIX_ARGS=()
if command -v brew >/dev/null 2>&1; then
  BREW_PREFIX="$(brew --prefix 2>/dev/null || true)"
  [[ -n "$BREW_PREFIX" ]] && PREFIX_ARGS+=("-DCMAKE_PREFIX_PATH=$BREW_PREFIX")
fi

# --------------------------------------------------------------------------
# Stale-cache / generator-mismatch recovery
# --------------------------------------------------------------------------
if [[ "$CLEAN" -eq 1 ]]; then
  info "Removing build directory: $BUILD_DIR"
  rm -rf "$ROOT/$BUILD_DIR"
elif [[ -f "$ROOT/$BUILD_DIR/CMakeCache.txt" ]]; then
  cached_gen="$(grep -E '^CMAKE_GENERATOR:' "$ROOT/$BUILD_DIR/CMakeCache.txt" | cut -d= -f2- || true)"
  cached_src="$(grep -E '^CMAKE_HOME_DIRECTORY:' "$ROOT/$BUILD_DIR/CMakeCache.txt" | cut -d= -f2- || true)"
  if [[ -n "$cached_gen" && "$cached_gen" != "$GENERATOR" ]]; then
    warn "Cached generator '$cached_gen' != '$GENERATOR' — wiping $BUILD_DIR to avoid a generator-mismatch error."
    rm -rf "$ROOT/$BUILD_DIR"
  elif [[ -n "$cached_src" && "$cached_src" != "$ROOT" ]]; then
    warn "Cached source dir '$cached_src' != '$ROOT' — wiping $BUILD_DIR."
    rm -rf "$ROOT/$BUILD_DIR"
  fi
fi

# --------------------------------------------------------------------------
# Configure
# --------------------------------------------------------------------------
configure_hints() {
  warn "CMake configuration failed. Common causes:"
  cat >&2 <<'EOF'
    • Boost not found      -> brew install boost   (or run ./setup.sh)
    • Python3 not found    -> brew install python  (find_package(Python3 REQUIRED))
    • Submodules missing   -> git submodule update --init --recursive
    • Stale cache          -> ./build.sh --clean
EOF
}

info "Configuring ($BUILD_TYPE, generator: $GENERATOR) -> $BUILD_DIR"
# CMAKE_POLICY_VERSION_MINIMUM=3.5 keeps CMake 4.x happy with vendored
# libraries that still declare cmake_minimum_required(VERSION < 3.5).
if ! cmake -S "$ROOT" -B "$ROOT/$BUILD_DIR" -G "$GENERATOR" \
      -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
      -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
      "${PREFIX_ARGS[@]}"; then
  configure_hints
  die "Configure step failed (see hints above)."
fi

# --------------------------------------------------------------------------
# Build
# --------------------------------------------------------------------------
BUILD_ARGS=(--build "$ROOT/$BUILD_DIR" --config "$BUILD_TYPE" --parallel "$JOBS")
if [[ -n "$TARGET" ]]; then
  BUILD_ARGS+=(--target "$TARGET")
  info "Building target '$TARGET' with $JOBS jobs"
else
  info "Building all targets with $JOBS jobs"
fi

if ! cmake "${BUILD_ARGS[@]}"; then
  warn "Build failed. If this is the first build after changing toolchains, try: ./build.sh --clean"
  die "Compilation failed (see errors above)."
fi

# --------------------------------------------------------------------------
# Report artefacts
# --------------------------------------------------------------------------
ART="$ROOT/$BUILD_DIR/Quotient_artefacts/$BUILD_TYPE"
printf "\n${C_GREEN}Build complete (%s).${C_RESET}\n" "$BUILD_TYPE"
if [[ -d "$ART" ]]; then
  echo "Artefacts:"
  for kind in VST3 AU Standalone; do
    if compgen -G "$ART/$kind/*" >/dev/null; then
      for p in "$ART/$kind/"*; do
        printf "    %-11s %s\n" "$kind" "$p"
      done
    fi
  done
fi
