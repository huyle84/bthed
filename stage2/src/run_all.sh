#!/bin/sh
# Stage 2a batch. Usage: ./run_all.sh   (after `make`)
set -e
R=../results; mkdir -p $R
./bthed_fhe2 20 8 5 3 1 > $R/fhe_mldsa44_L8_n5t3.csv
./bthed_fhe2 10 8 7 4 2 > $R/fhe_mldsa44_L8_n7t4.csv
./bthed_fhe3 10 8 5 3 3 > $R/fhe_mldsa65_L8_n5t3.csv
./bthed_fhe5 10 8 5 3 4 > $R/fhe_mldsa87_L8_n5t3.csv
# Stage 2b: one comparison (inputs: multiples of 2^-alpha), sign ~ f3^df o g3^dg
./bthed_cmp 8 4 2 1 > $R/cmp_alpha8.csv
./bthed_cmp 10 5 2 1 > $R/cmp_alpha10.csv
./bthed_cmp 12 6 2 1 > $R/cmp_alpha12.csv
{ uname -srm; grep -m1 'model name' /proc/cpuinfo 2>/dev/null; echo "cores: $(nproc)"; ${CC:-cc} --version | head -1; echo "flags: -O3 -march=native"; } > $R/environment.txt
