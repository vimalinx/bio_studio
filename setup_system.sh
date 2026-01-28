#!/bin/bash
set -e

echo "🔧 Bio Studio System Setup (Arch Linux)"
echo "======================================"

# 1. System Tools
echo "📦 [1/3] Installing system tools via Pacman..."
echo "   (You may be asked for your sudo password)"
sudo pacman -S --needed --noconfirm \
    base-devel cmake \
    tree jq ripgrep fd unzip zip \
    net-tools htop git wget curl

# 2. Bioinformatics Infrastructure
echo "🧬 [2/4] Installing Bioinformatics System Tools..."
sudo pacman -S --needed --noconfirm \
    apptainer \
    graphviz \
    gnuplot \
    imagemagick \
    openmpi \
    r \
    python-numpy

# 3. Fix/Install SDKMAN
echo "☕ [3/4] Setting up SDKMAN..."
export SDKMAN_DIR="$HOME/.sdkman"

# Check for broken install (folder exists but no bin)
if [ -d "$SDKMAN_DIR" ] && [ ! -f "$SDKMAN_DIR/bin/sdkman-init.sh" ]; then
    echo "   ⚠️ Detected incomplete SDKMAN installation. Re-installing..."
    rm -rf "$SDKMAN_DIR"
fi

if [ ! -d "$SDKMAN_DIR" ]; then
    curl -s "https://get.sdkman.io" | bash
else
    echo "   ✅ SDKMAN is intact."
fi

# Try to load SDKMAN to install Java
# We need to source it dynamically
if [ -f "$SDKMAN_DIR/bin/sdkman-init.sh" ]; then
    source "$SDKMAN_DIR/bin/sdkman-init.sh"
    
    echo "☕ [4/4] Installing Java Versions..."
    
    # Check if Java 17 is already installed
    if ! sdk list java | grep -q "17.0.10-tem.*installed"; then
        sdk install java 17.0.10-tem
    else
        echo "   ✅ Java 17 already installed"
    fi

    # Check if Java 8 is already installed
    if ! sdk list java | grep -q "8.0.402-tem.*installed"; then
        sdk install java 8.0.402-tem
    else
        echo "   ✅ Java 8 already installed"
    fi

    sdk default java 17.0.10-tem
else
    echo "❌ Error: SDKMAN failed to initialize. Please run 'curl -s https://get.sdkman.io | bash' manually."
    exit 1
fi

echo ""
echo "🎉 Setup Complete!"
echo "----------------------------------------"
echo "✅ System Tools: cmake, tree, ripgrep..."
echo "✅ Java Environment: 17 (Default), 8"
echo ""
echo "👉 IMPORTANT: Close and Re-open your terminal to make 'sdk' command available."
