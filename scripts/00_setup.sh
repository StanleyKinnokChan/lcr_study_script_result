#!/bin/bash
# Install fLPS 2.0 and Python dependencies
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

echo "=== Installing fLPS 2.0 ==="
cd "$PROJECT_DIR"

if [ ! -d "flps2" ]; then
    git clone https://github.com/pmharrison/flps2.git
fi

cd flps2

if [ ! -f "flps" ]; then
    # Extract tar.gz if not already done
    if [ ! -d "fLPS2programs" ]; then
        TARBALL=$(ls fLPS2programs*.tar.gz 2>/dev/null | head -1)
        if [ -z "$TARBALL" ]; then
            echo "ERROR: fLPS2programs.tar.gz not found in $(pwd)"
            exit 1
        fi
        echo "Extracting $TARBALL..."
        tar -xzf "$TARBALL"
    fi

    # Pick the right pre-built binary for this OS (no compilation needed)
    OS="$(uname -s)"
    case "$OS" in
        Darwin) BIN_DIR="fLPS2programs/bin/macosx" ;;
        Linux)  BIN_DIR="fLPS2programs/bin/linux"  ;;
        *)
            echo "ERROR: Unsupported OS: $OS"
            exit 1
            ;;
    esac

    SRC_BIN="$BIN_DIR/fLPS2"
    if [ ! -f "$SRC_BIN" ]; then
        echo "ERROR: Pre-built binary not found at $SRC_BIN"
        exit 1
    fi

    cp "$SRC_BIN" flps
    chmod +x flps
    echo "Copied $SRC_BIN -> flps2/flps"
fi

echo "fLPS 2.0 binary: $(pwd)/flps"
cd "$PROJECT_DIR"

echo ""
echo "=== Installing Python dependencies ==="
pip install requests pandas scipy matplotlib seaborn biopython --quiet

echo ""
echo "=== Setup complete ==="
echo "fLPS binary: $PROJECT_DIR/flps2/flps"
echo "Run scripts in order: 01 -> 02 -> 03 -> 04"
