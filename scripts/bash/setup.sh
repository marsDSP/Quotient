#!/usr/bin/env bash
#
# setup.sh — one-time environment setup for the Quotient JUCE plugin.
#
# This script makes a freshly-cloned checkout buildable on a new machine:
#   1. Verifies the core toolchain (git, a C/C++ compiler / Xcode CLT).
#   2. Installs build dependencies (cmake, ninja, boost, python, doxygen).
#   3. Initialises every git submodule under libs/.
#   4. Makes the helper scripts executable.
#
# After running this once, use ./scripts/bash/build.sh to configure and compile.
#
# Usage:
#   ./scripts/bash/setup.sh [--skip-deps] [--skip-submodules] [-h|--help]
#
set -euo pipefail

# --------------------------------------------------------------------------
# Paths / constants
# --------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# This script lives in <repo>/scripts/bash, so the project root (where
# CMakeLists.txt, libs/ and .gitmodules live) is two directories up.
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$ROOT"

# Public upstream JUCE, used only as a fallback if submodule init cannot
# populate libs/JUCE (see the project rule on JUCE setup).
JUCE_FALLBACK_URL="https://github.com/juce-framework/JUCE.git"

BREW_PACKAGES=(cmake ninja boost doxygen pkg-config)

SKIP_DEPS=0
SKIP_SUBMODULES=0

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
setup.sh — one-time environment setup for the Quotient JUCE plugin.

Prepares a freshly-cloned checkout to build:
  1. Verifies the core toolchain (git, Xcode Command Line Tools on macOS).
  2. Installs build dependencies (cmake, ninja, boost, python, doxygen).
  3. Initialises every git submodule under libs/ (JUCE, Catch2, kfr,
     pluginval, gcem, tracy, simde, xsimd).
  4. Makes the helper scripts executable.

After running this once, use ./scripts/bash/build.sh to configure and compile.

Usage:
  ./scripts/bash/setup.sh [--skip-deps] [--skip-submodules] [-h|--help]

Options:
  --skip-deps          Do not install/verify build dependencies.
  --skip-submodules    Do not initialise git submodules.
  -h, --help           Show this help.
EOF
  exit 0
}

# --------------------------------------------------------------------------
# Argument parsing
# --------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --skip-deps)       SKIP_DEPS=1 ;;
    --skip-submodules) SKIP_SUBMODULES=1 ;;
    -h|--help)         usage ;;
    *) die "Unknown option: $1 (try --help)" ;;
  esac
  shift
done

OS="$(uname -s)"

# --------------------------------------------------------------------------
# 1. Core toolchain
# --------------------------------------------------------------------------
info "Checking core toolchain"
command -v git >/dev/null 2>&1 || die "git is required but was not found on PATH."
ok "git: $(git --version | awk '{print $3}')"

if [[ "$OS" == "Darwin" ]]; then
  if ! xcode-select -p >/dev/null 2>&1; then
    warn "Xcode Command Line Tools are not installed."
    info "Launching the installer (a GUI dialog will appear)…"
    xcode-select --install || true
    die "Re-run ./setup.sh once the Command Line Tools finish installing."
  fi
  ok "Xcode command line tools: $(xcode-select -p)"
fi

# --------------------------------------------------------------------------
# 2. Build dependencies
# --------------------------------------------------------------------------
ensure_python3() {
  if command -v python3 >/dev/null 2>&1; then
    ok "python3: $(python3 --version 2>&1 | awk '{print $2}')"
  else
    warn "python3 not found — CMakeLists.txt requires it (find_package(Python3 REQUIRED))."
    return 1
  fi
}

if [[ "$SKIP_DEPS" -eq 1 ]]; then
  info "Skipping dependency installation (--skip-deps)"
elif [[ "$OS" == "Darwin" ]]; then
  info "Installing build dependencies via Homebrew"
  if ! command -v brew >/dev/null 2>&1; then
    warn "Homebrew is not installed."
    echo '    Install it with:'
    echo '      /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"'
    die "Re-run ./setup.sh after installing Homebrew."
  fi
  for pkg in "${BREW_PACKAGES[@]}"; do
    if brew list --versions "$pkg" >/dev/null 2>&1; then
      ok "$pkg already installed"
    else
      info "brew install $pkg"
      brew install "$pkg"
    fi
  done
  # python3: prefer an existing interpreter, otherwise install the brew formula.
  ensure_python3 || { info "brew install python"; brew install python; ensure_python3 || die "python3 still not found after install."; }
elif [[ "$OS" == "Linux" ]]; then
  info "Linux detected — install the following with your package manager, then re-run:"
  if   command -v apt-get >/dev/null 2>&1; then
    echo "    sudo apt-get update && sudo apt-get install -y build-essential cmake ninja-build libboost-all-dev python3 python3-dev doxygen pkg-config libasound2-dev libx11-dev libxext-dev libxrandr-dev libxinerama-dev libxcursor-dev libfreetype6-dev libfontconfig1-dev libgl1-mesa-dev libcurl4-openssl-dev"
  elif command -v dnf >/dev/null 2>&1; then
    echo "    sudo dnf install -y gcc-c++ cmake ninja-build boost-devel python3 python3-devel doxygen pkgconf alsa-lib-devel libX11-devel libXext-devel libXrandr-devel libXinerama-devel libXcursor-devel freetype-devel fontconfig-devel mesa-libGL-devel libcurl-devel"
  elif command -v pacman >/dev/null 2>&1; then
    echo "    sudo pacman -S --needed base-devel cmake ninja boost python doxygen pkgconf alsa-lib libx11 libxext libxrandr libxinerama libxcursor freetype2 fontconfig mesa curl"
  else
    echo "    (cmake ninja boost python3 doxygen pkg-config + JUCE Linux dev headers)"
  fi
  ensure_python3 || warn "Install python3 before building."
else
  warn "Unrecognised OS '$OS' — install cmake, ninja, boost and python3 manually."
fi

# --------------------------------------------------------------------------
# 3. Git submodules
# --------------------------------------------------------------------------
if [[ "$SKIP_SUBMODULES" -eq 1 ]]; then
  info "Skipping submodule initialisation (--skip-submodules)"
elif [[ -f "$ROOT/.gitmodules" ]]; then
  info "Initialising git submodules under libs/"
  git submodule sync --recursive
  if ! git submodule update --init --recursive; then
    warn "git submodule update failed — see messages above (private fork or network issue?)."
  fi

  # The build hard-requires JUCE (add_subdirectory(libs/JUCE)). If it is still
  # missing, fall back to cloning upstream JUCE so the repo can still build.
  if [[ ! -f "$ROOT/libs/JUCE/CMakeLists.txt" ]]; then
    warn "libs/JUCE is empty after submodule init."
    info "Falling back to cloning JUCE from $JUCE_FALLBACK_URL"
    rm -rf "$ROOT/libs/JUCE"
    git clone --depth 1 "$JUCE_FALLBACK_URL" "$ROOT/libs/JUCE" \
      || die "Failed to clone JUCE into libs/JUCE."
  fi
  ok "Submodules ready"
else
  warn ".gitmodules not found — skipping submodule init."
fi

# --------------------------------------------------------------------------
# 4. Helper scripts: ensure they are executable + optional python venv
# --------------------------------------------------------------------------
info "Configuring helper scripts"
if compgen -G "$ROOT/scripts/bash/*.sh" >/dev/null; then
  chmod +x "$ROOT"/scripts/bash/*.sh
  ok "scripts/bash/*.sh marked executable"
fi

# If a Python requirements file ships with the repo, set up a local venv.
for req in "$ROOT/requirements.txt" "$ROOT/scripts/python/requirements.txt"; do
  if [[ -f "$req" ]]; then
    info "Creating Python virtualenv for $req"
    python3 -m venv "$ROOT/.venv"
    # shellcheck disable=SC1091
    source "$ROOT/.venv/bin/activate"
    pip install --quiet --upgrade pip
    pip install --quiet -r "$req"
    deactivate
    ok "Python venv ready at .venv (activate with: source .venv/bin/activate)"
    break
  fi
done

# --------------------------------------------------------------------------
# Done
# --------------------------------------------------------------------------
printf "\n${C_GREEN}Setup complete.${C_RESET} Next:\n"
printf "    ./scripts/bash/build.sh            # debug build (default)\n"
printf "    ./scripts/bash/build.sh --release  # optimised build\n"
printf "    ./scripts/bash/build.sh --help     # all options\n"
