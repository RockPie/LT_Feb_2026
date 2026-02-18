#!/bin/bash
# Clean build script that ensures conda is not active during compilation

# Save current directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "===== Checking conda status ====="
if [ -n "$CONDA_PREFIX" ]; then
    echo "ERROR: Conda is still active!"
    echo "Please run 'conda deactivate' until conda is completely deactivated,"
    echo "then run this script again."
    exit 1
fi

echo "===== Conda is deactivated ====="

# Set up ROOT environment
export ROOTSYS=/home/shihai/sw/root/root_install
export PATH=$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH

# Ensure we're using system compilers, not conda compilers
export CC=/usr/bin/gcc
export CXX=/usr/bin/g++

echo "===== Cleaning old build ====="
rm -rf build
mkdir build
cd build

echo "===== Running CMake ====="
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_COMPILER=/usr/bin/gcc \
      -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
      ..

if [ $? -ne 0 ]; then
    echo "ERROR: CMake configuration failed"
    exit 1
fi

echo "===== Building ====="
make -j$(nproc)

if [ $? -ne 0 ]; then
    echo "ERROR: Build failed"
    exit 1
fi

echo "===== Build successful ====="
echo "You can now run executables directly without the wrapper script:"
echo "  ./build/103_Rootifier_10G -f data/Run300.h2g -o dump/103_Rootifier_10G/Run300_rootified.root"
