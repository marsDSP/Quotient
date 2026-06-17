#!/bin/bash
# copy.sh — POST_BUILD helper that copies the freshly-built VST3 bundle
# into a local "test" folder for quick scanning by a DAW/host.
BUNDLE_PATH="$1"
if [ -z "$BUNDLE_PATH" ]; then
    BUNDLE_PATH="$HOME/CLionProjects/Quotient/cmake-build-debug/Quotient_artefacts/Debug/VST3/Quotient.vst3"
fi

DEST_DIR="$HOME/Desktop/vst test"

if [ ! -d "$BUNDLE_PATH" ]; then
    echo "copy_vst3.sh: VST3 bundle not found at '$BUNDLE_PATH'; skipping copy."
    exit 0
fi

VST3_NAME="$(basename "$BUNDLE_PATH")"

# Create destination directory if it doesn't exist
if ! mkdir -p "$DEST_DIR"; then
    echo "copy_vst3.sh: could not create '$DEST_DIR'; skipping copy."
    exit 0
fi

# Copy the .vst3 bundle (it's a directory on macOS).
# shellcheck disable=SC2115
rm -rf "$DEST_DIR/$VST3_NAME" # Remove old version first to ensure clean copy
if cp -R "$BUNDLE_PATH" "$DEST_DIR/"; then
    echo "Successfully copied $VST3_NAME to $DEST_DIR"
else
    echo "copy_vst3.sh: failed to copy $VST3_NAME (non-fatal)."
fi
exit 0
