#!/usr/bin/env bash
# The system libraries GPUI builds against on Debian and Ubuntu.
set -euo pipefail
apt-get update -q
apt-get install -qy --no-install-recommends \
    build-essential clang pkg-config cmake \
    libxkbcommon-dev libxkbcommon-x11-dev \
    libwayland-dev libx11-dev libx11-xcb-dev libxcb1-dev libxi-dev \
    libfontconfig1-dev libfreetype6-dev \
    libvulkan-dev libvulkan1 mesa-vulkan-drivers \
    libasound2-dev libssl-dev libzstd-dev \
    imagemagick
