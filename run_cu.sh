#!/bin/bash

# Generate input files
cd ~/seed-extension/local_align
gcc -o rands rand_str.cpp
./rands

# Compile
cd ~/seed-extension/local_align
nvcc -o swout_cu swalign.cu

# Run
./swout_cu
