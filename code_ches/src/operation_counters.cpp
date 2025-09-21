#include <stdio.h>

#include "operation_counters.h"

operation_counters_t FHEOperationCounters = {0};

void operation_counters_init() {
    FHEOperationCounters.counter__gsw_automorphism_RNSc_half_gsw = 0;
    FHEOperationCounters.counter__gsw_RNS_to_half_RNSc_gsw = 0;
    FHEOperationCounters.counter__gsw_mul_RNSc_half_gsw = 0;
    FHEOperationCounters.counter__rlwe_automorphism_RNSc = 0;
    FHEOperationCounters.counter__gsw_mul_RNSc_rlwe = 0;
}

void operation_counters_increment_gsw_automorphism_RNSc_half_gsw() {
    #pragma omp atomic
    FHEOperationCounters.counter__gsw_automorphism_RNSc_half_gsw++;
}

void operation_counters_increment_gsw_RNS_to_half_RNSc_gsw() {
    #pragma omp atomic
    FHEOperationCounters.counter__gsw_RNS_to_half_RNSc_gsw++;
    // printf("[operation_counters_increment_gsw_mul_RNSc_half_gsw: %ld]\n", FHEOperationCounters.counter__gsw_RNS_to_half_RNSc_gsw);
}

void operation_counters_increment_gsw_mul_RNSc_half_gsw() {
    #pragma omp atomic
    FHEOperationCounters.counter__gsw_mul_RNSc_half_gsw++;
}

void operation_counters_increment_rlwe_automorphism_RNSc() {
    #pragma omp atomic
    FHEOperationCounters.counter__rlwe_automorphism_RNSc++;
}

void operation_counters_increment_gsw_mul_RNSc_rlwe() {
    #pragma omp atomic
    FHEOperationCounters.counter__gsw_mul_RNSc_rlwe++;
}

void operation_counters_show() {
    fprintf(stderr, "FHE operation counters ============\n");
    fprintf(stderr, "    gsw_automorphism_RNSc_half_gsw: ");
    fprintf(stderr, "%ld\n", FHEOperationCounters.counter__gsw_automorphism_RNSc_half_gsw);
    fprintf(stderr, "    gsw_RNS_to_half_RNSc_gsw: ");
    fprintf(stderr, "%ld\n", FHEOperationCounters.counter__gsw_RNS_to_half_RNSc_gsw);
    fprintf(stderr, "    gsw_mul_RNSc_half_gsw: ");
    fprintf(stderr, "%ld\n", FHEOperationCounters.counter__gsw_mul_RNSc_half_gsw);
    fprintf(stderr, "    rlwe_automorphism_RNSc: ");
    fprintf(stderr, "%ld\n", FHEOperationCounters.counter__rlwe_automorphism_RNSc);
    fprintf(stderr, "    gsw_mul_RNSc_rlwe: ");
    fprintf(stderr, "%ld\n", FHEOperationCounters.counter__gsw_mul_RNSc_rlwe);
    fprintf(stderr, "===================================\n");

}