#!/bin/bash

if [[ -z $1 ]];
then
    echo "[-] Usage $0 <output_csv_file>"
    exit 1
fi

set -e # Stop on any error
set -x

output_csv_file=$1

export OMP_NUM_THREADS=1 # Used because we want accurate performance results
ulimit -s 1000000  # Needed for N >= 2048

make hexl
make

N_BATCHES=1

./main_amortized print_header | tee -a $output_csv_file

for i in `seq 1 $N_BATCHES`; do
    echo "---- Running batch $i -----"

    ./main_amortized 12289 8 1024 4 6 3 24 256 | tee -a $output_csv_file
    ./main_amortized 7681 8 1024 4 6 3 24 256 | tee -a $output_csv_file
    ./main_amortized 7937 8 1024 4 6 3 24 256 | tee -a $output_csv_file
    ./main_amortized 16001 8 1024 4 6 3 24 256 | tee -a $output_csv_file
    ./main_amortized 7681 8 2048 4 6 3 27 52 | tee -a $output_csv_file
    ./main_amortized 7937 8 2048 4 6 3 27 52 | tee -a $output_csv_file
    ./main_amortized 12289 8 2048 4 6 3 27 52 | tee -a $output_csv_file
    ./main_amortized 15361  8 2048 4 6 3 27 52 | tee -a $output_csv_file

done

set +x
