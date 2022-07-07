#!/bin/bash -e

opt=" -w -Wextra -pedantic -g -O2 -fallow-argument-mismatch"

cmake \
-DCMAKE_C_COMPILER=mpicc \
-DCMAKE_CXX_COMPILER=mpicxx \
-DCMAKE_Fortran_COMPILER=gfortran \
-DCMAKE_Fortran_FLAGS="-fallow-argument-mismatch" \
-DCMAKE_C_FLAGS="${opt}" \
-DCMAKE_CXX_FLAGS="${opt}" \
-DCMAKE_BUILD_TYPE=Debug \
-DPHASTA_INCOMPRESSIBLE=OFF \
-DPHASTA_COMPRESSIBLE=ON \
..

make -j8
