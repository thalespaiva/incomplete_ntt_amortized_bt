#pragma once

#include <stddef.h>
#include <stdint.h>

typedef struct operation_counters_s {
    size_t counter__gsw_automorphism_RNSc_half_gsw;
    size_t counter__gsw_RNS_to_half_RNSc_gsw;
    size_t counter__gsw_mul_RNSc_half_gsw;
    size_t counter__rlwe_automorphism_RNSc;
    size_t counter__gsw_mul_RNSc_rlwe;
} operation_counters_t;

extern operation_counters_t FHEOperationCounters;

void operation_counters_init();

void operation_counters_increment_gsw_automorphism_RNSc_half_gsw(); // OK
void operation_counters_increment_gsw_RNS_to_half_RNSc_gsw();       // OK
void operation_counters_increment_gsw_mul_RNSc_half_gsw();          // OK
void operation_counters_increment_rlwe_automorphism_RNSc();         // OK
void operation_counters_increment_gsw_mul_RNSc_rlwe();              // OK
void operation_counters_show();
