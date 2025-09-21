#include <climits>

#include "rns-gsw.h"

uint64_t MODULUS_P = 0;
uint64_t *discrete_log_base_gen_mod_q = NULL;
uint64_t Z_P_STAR_GENERATOR = 0;

void set_modulus_p(uint64_t modulus_p) {
    assert(intel::hexl::IsPrime(modulus_p));

    MODULUS_P = modulus_p;
    fprintf(stderr, "Setting MODULUS_P as %ld\n", MODULUS_P);

    discrete_log_base_gen_mod_q = (uint64_t *) malloc((modulus_p - 1) * sizeof(*discrete_log_base_gen_mod_q));

    Z_P_STAR_GENERATOR = intel::hexl::GeneratePrimitiveRoot(modulus_p - 1, modulus_p);
    for (size_t i = 0; i < modulus_p - 1; i++) {
        uint64_t power = intel::hexl::PowMod(Z_P_STAR_GENERATOR, i, modulus_p);

        discrete_log_base_gen_mod_q[power - 1] = i;
    }
}
