#!/bin/bash

set -e # Stop on any error

# Some lists of primes that can be useful
# primes12_bits_for_N1024="2689, 3329, 3457"
# primes13_bits_for_N1024="4481, 4993, 6529, 7297, 7681, 7937"
# primes14_bits_for_N1024="9473, 9601, 9857, 10369, 10753, 11393, 11777, 12161, 12289, 13313, 13441, 13697, 14081, 14593, 15233, 15361, 16001"


# Primes that will be used in this experiment
selected_primes_for_N1024="7297, 7681, 7937, 9473, 10753, 12289, 14081, 16001"

make clean
make main_amortized

# We need to increase the stack otherwise we get a seg fault
ulimit -s 10000

# Only print the header to the output
./main_amortized print_header

for MODULUS_P in $selected_primes_for_N1024; do
    for i in `seq 1 10`; do
        # "Usage: %s <prime_p> <n_message_bits> <N> <NTT_INCOMPLETENESS> <INTT_QUADRATIC_LEVEL> <l> <p_star_size> <key_hw>"
        ./main_amortized $MODULUS_P 11 1024 4 6 3 24 256
    done
done
