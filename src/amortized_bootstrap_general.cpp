Copyright 2023 Antonio Guimarães et al.
Copyright 2025 LG Electronics, Inc.
SPDX-License-Identifier: Apache-2.0

#include "rns-gsw.h"
#include "benchmark_util.h"
#include "perf_counters.h"
#include "omp.h"

#define INTT_MAX_DEPTH 2
#define INTT_M 1

#define SHRINK_LEVELS 0

// #define MEASURE_CACHE_MISSES
#ifdef MEASURE_CACHE_MISSES
#define perf_counters_start_if_measure_cache_misses_is_set() perf_counters_start();
#define perf_counters_stop_if_measure_cache_misses_is_set() perf_counters_stop();
#else
#define perf_counters_start_if_measure_cache_misses_is_set() {};
#define perf_counters_stop_if_measure_cache_misses_is_set() {};
#endif

// Note: Decomposed automorphism will be added in a future version
#ifdef USE_DECOMP_AUTOMORPHISM
#define MACRO_rlwe_automorphism(OUT, IN, GEN) rlwe_decomp_automorphism_RNSc(OUT, IN, GEN, aut_ksk);
#define MACRO_hgsw_automorphism(OUT, IN, GEN) gsw_decomp_automorphism_RNSc_half_gsw(OUT, IN, GEN, aut_ksk);
#else
#define MACRO_rlwe_automorphism(OUT, IN, GEN) rlwe_automorphism_RNSc(OUT, IN, GEN, aut_ksk[GEN]);
#define MACRO_hgsw_automorphism(OUT, IN, GEN) gsw_automorphism_RNSc_half_gsw(OUT, IN, GEN, aut_ksk[GEN]);
#endif

#define MAX_THREADS 128

extern int INTT_QUADRATIC_LEVEL;
extern int NTT_INCOMPLETENESS;

static RNSc_half_GSW * acc, tmp, tmp_for_threads[MAX_THREADS];

// The intt algorithm does not run in place, so we need a buffer for temporary results.
void setup_tmp_buffer(uint64_t n, uint64_t Q, uint64_t l, intel::hexl::NTT ** ntt){
  acc = (RNSc_half_GSW *) gsw_alloc_half_RNS_sample_array(n, Q, l, ntt);
  tmp = (RNSc_half_GSW) gsw_alloc_half_RNS_sample(Q, l, ntt);

  // TODO: this is not the ideal place to put this
  if (omp_get_max_threads() > MAX_THREADS) {
    omp_set_num_threads(MAX_THREADS);
  }

  for (int i = 0; i < MAX_THREADS; i++) {
    tmp_for_threads[i] = (RNSc_half_GSW) gsw_alloc_half_RNS_sample(Q, l, ntt);
  }
}


// Implements the dot product with GSW accumulators
void rns_gsw_dot_product(RNSc_half_GSW out_acc,
                         RNS_GSW * encrypted_coefficients,
                         uint64_t plaintext_coefficients[],
                         size_t length,
                         RNS_RLWE_KS_Key * aut_ksk,
                         RNS_RLWE_KS_Key * priv_ksk) {

  const uint64_t Q = encrypted_coefficients[0]->samples[0]->b->Q;

  size_t nonzero_indexes[length + 1] = {0};
  size_t n_nonzero_indexes = 0;
  for (size_t i = 0; i < length; i++) {
    if (plaintext_coefficients[i] != 0) {
      nonzero_indexes[n_nonzero_indexes++] = i;
    }
  }
  // We need to set `nonzero_indexes[n_nonzero_indexes]` in a way that is compatible to Line 1 (Algorithm 7 EvalScalarProd)
  // This ensures that `plaintext_inverse_coefficients[nonzero_indexes[n_nonzero_indexes]] = 1`.
  nonzero_indexes[n_nonzero_indexes] = length;

  // if (n_nonzero_indexes != length) {
  //   fprintf(stderr, "[* rns_gsw_dot_product] ");
  //   fprintf(stderr, "Found %ld non-null coefficients\n", n_nonzero_indexes);
  // }
  assert(n_nonzero_indexes >= 1);


  uint64_t plaintext_inverse_coefficients[length + 1] = {0};

  for (size_t i_nz = 0; i_nz < n_nonzero_indexes; i_nz++) {
    size_t i = nonzero_indexes[i_nz];
    plaintext_inverse_coefficients[i] = intel::hexl::InverseMod(plaintext_coefficients[i], Q);
  }
  plaintext_inverse_coefficients[length] = 1;  // From Line 1 (Algorithm 7 EvalScalarProd)


  gsw_RNS_to_half_RNSc_gsw(out_acc, encrypted_coefficients[nonzero_indexes[0]]);
  const uint64_t gen = intel::hexl::MultiplyMod(plaintext_coefficients[nonzero_indexes[0]],
                                                plaintext_inverse_coefficients[nonzero_indexes[1]], Q);


  MACRO_hgsw_automorphism(out_acc, out_acc, gen);

  for (size_t j = 1; j < n_nonzero_indexes; j++) {
    // In the last iteration `j + 1 = n_nonzero_indexes`, but we remember we already ensured that:
    // plaintext_inverse_coefficients[nonzero_indexes[n_nonzero_indexes]] = 1;
    const uint64_t gen = intel::hexl::MultiplyMod(plaintext_coefficients[nonzero_indexes[j]],
                                                  plaintext_inverse_coefficients[nonzero_indexes[j + 1]], Q);
    gsw_mul_RNSc_half_gsw(tmp_for_threads[omp_get_thread_num()], encrypted_coefficients[nonzero_indexes[j]], out_acc);

    MACRO_hgsw_automorphism(out_acc, tmp_for_threads[omp_get_thread_num()], gen);
  }
}


// Implements the dot product with RLWE accumulators.
void rnsc_rlwe_dot_product(RNSc_RLWE out_tv,
                           RNS_GSW * encrypted_coefficients,
                           uint64_t plaintext_coefficients[],
                           size_t length,
                           RNS_RLWE_KS_Key * aut_ksk,
                           RNS_RLWE_KS_Key * priv_ksk) {

  const uint64_t Q = encrypted_coefficients[0]->samples[0]->b->Q;

  size_t nonzero_indexes[length + 1] = {0};
  size_t n_nonzero_indexes = 0;
  for (size_t i = 0; i < length; i++) {
    if (plaintext_coefficients[i] != 0) {
      nonzero_indexes[n_nonzero_indexes++] = i;
    }
  }
  // We need to set `nonzero_indexes[n_nonzero_indexes]` in a way that is compatible to Line 1 (Algorithm 7 EvalScalarProd)
  // This ensures that `plaintext_inverse_coefficients[nonzero_indexes[n_nonzero_indexes]] = 1`.
  nonzero_indexes[n_nonzero_indexes] = length;

  // if (n_nonzero_indexes != length) {
  //   fprintf(stderr, "[* rnsc_rlwe_dot_product] ");
  //   fprintf(stderr, "Found %ld non-null coefficients\n", n_nonzero_indexes);
  // }
  assert(n_nonzero_indexes >= 1);


  uint64_t plaintext_inverse_coefficients[length + 1] = {0};

  for (size_t i_nz = 0; i_nz < n_nonzero_indexes; i_nz++) {
    size_t i = nonzero_indexes[i_nz];
    plaintext_inverse_coefficients[i] = intel::hexl::InverseMod(plaintext_coefficients[i], Q);
  }
  plaintext_inverse_coefficients[length] = 1;  // From Line 1 (Algorithm 7 EvalScalarProd)

  // The dot product for the RLWE case is different -- see Section 3.7.
  const uint64_t gen = intel::hexl::InverseMod(plaintext_coefficients[nonzero_indexes[0]], Q);
  MACRO_rlwe_automorphism(out_tv, out_tv, gen);

  for (size_t j = 0; j < n_nonzero_indexes; j++) {
    // In the last iteration `j + 1 = n_nonzero_indexes`, but we remember we already ensured that:
    // plaintext_inverse_coefficients[nonzero_indexes[n_nonzero_indexes]] = 1;
    const uint64_t gen = intel::hexl::MultiplyMod(plaintext_coefficients[nonzero_indexes[j]],
                                                  plaintext_inverse_coefficients[nonzero_indexes[j + 1]], Q);

    gsw_mul_RNSc_rlwe(tmp_for_threads[omp_get_thread_num()]->samples[0], encrypted_coefficients[nonzero_indexes[j]], out_tv);
    MACRO_rlwe_automorphism(out_tv, tmp_for_threads[omp_get_thread_num()]->samples[0], gen);
  }
}


void get_middle_intt_matrix_column(uint64_t column[], int idx, uint64_t N, uint64_t p, uint64_t rho, uint64_t ntt_incompleteness) {

  uint64_t rou = intel::hexl::MinimalPrimitiveRoot((2*N) >> ntt_incompleteness, p);
  size_t rou_order = (2*N) >> ntt_incompleteness;

  size_t max_n_levels = ilog2(N);

  uint64_t one_hot_vector[N] = {0};
  one_hot_vector[idx] = 1;

  uint64_t outA[N] = {0};
  nc_ntt_inverse_leveled(outA, one_hot_vector, rou, rou_order, p, N, ntt_incompleteness);

  nc_ntt_forward_leveled(column, outA, rou, rou_order, p, N, max_n_levels - rho);
}

void get_part2_intt_matrix_column(uint64_t column[], int idx, uint64_t N, uint64_t p, uint64_t rho, uint64_t ntt_incompleteness) {

  uint64_t rou = intel::hexl::MinimalPrimitiveRoot((2*N) >> ntt_incompleteness, p);
  size_t rou_order = (2*N) >> ntt_incompleteness;

  size_t max_n_levels = ilog2(N);

  uint64_t one_hot_vector[N] = {0};
  one_hot_vector[idx] = 1;

  nc_ntt_inverse_leveled(column, one_hot_vector, rou, rou_order, p, N, max_n_levels - rho);
}

// Returns the twiddle factors for the incomplete NTT for some incompleteness
// level that are used to define the rings over for the base multiplications.
// The length of the output is (n >> incompleteness)
void get_twiddles_for_basemul(uint64_t *out, uint64_t rou_2nth,
                              uint64_t Q, uint64_t n, int incompleteness) {

  uint64_t rou = intel::hexl::MinimalPrimitiveRoot((2*n) >> incompleteness, Q);

  uint64_t ws[n];
  calcRootsOfUnit(ws, rou, Q, n);

  int n_levels = ilog2(n) - incompleteness;
  int k = ilog2(n);

  int i = n_levels - 1;
  int d = 1 << (k - 1 - i);

  int idx_twiddle = (n >> incompleteness) - 1;
  for (int j = 0; j < (1 << i); j++) {
    int phi_idx = bit_rev[bit_rev_index_for_n(n >> incompleteness)][(1 << i) + j];
    uint64_t phi = intel::hexl::InverseMod(ws[phi_idx], Q);
    out[idx_twiddle--] = phi;
    out[idx_twiddle--] = intel::hexl::SubUIntMod(0, phi, Q);
  }

}

// Implements the computation of Part 1 of the INTT.
void intt_negacyclic_part1_quadratic(RNS_GSW * p, uint64_t * in_consts, RNSc_half_GSW * acc,
                                     uint64_t * w, uint64_t n,
                                     uint64_t N_main,
                                     RNS_RLWE_KS_Key * aut_ksk, RNS_RLWE_KS_Key * priv_ksk) {
  const uint64_t Q = p[0]->samples[0]->b->Q;

  // Build the incomplete INTT matrix corresponding to the part 1 (quadratic) step of the INTT computation
  uint64_t trivial_ntt_matrix[N_main][N_main] = {0};

  MEASURE_TIME("", 1, "Build-only Matrix 1", {
    for (size_t i = 0; i < N_main; i++) {
      uint64_t column_i[N_main];
      get_middle_intt_matrix_column(column_i, i, N_main, Q, INTT_QUADRATIC_LEVEL, NTT_INCOMPLETENESS);
      for (size_t j = 0; j < N_main; j++)
        trivial_ntt_matrix[j][i] = column_i[j];
    }
  })

  // Compute the `n >> NTT_INCOMPLETENESS` twiddles for the basemul operation
  uint64_t basemul_twiddles[n >> NTT_INCOMPLETENESS];
  fprintf(stderr, "w[1] = %ld\n", w[1]);
  get_twiddles_for_basemul(basemul_twiddles, w[1], Q, N_main, NTT_INCOMPLETENESS);
  fprintf(stderr, "twiddles for basemul: ");
  for (size_t i = 0; i < (n >> NTT_INCOMPLETENESS); i++) {
    fprintf(stderr, "%ld, ", basemul_twiddles[i]);
  }
  fprintf(stderr, "\n");

  fflush(stderr);

  #pragma omp parallel for
  for (size_t i = 0; i < n; i++) {
    uint64_t plaintext_coefficients[n] = {0};

    size_t basemul_dim = 1 << NTT_INCOMPLETENESS;

    for (size_t j = 0; j < n; j += basemul_dim) {
      uint64_t omega = basemul_twiddles[j / basemul_dim];  // Remember that `j` is always even

      for (size_t col = 0; col < basemul_dim; col++) {
        plaintext_coefficients[j + col] = 0;

        for (size_t row = 0; row < basemul_dim; row++) {
          size_t idx = (((row - col) % basemul_dim) + basemul_dim) % basemul_dim;

          if (row >= col) {
            plaintext_coefficients[j + col] += intel::hexl::MultiplyMod(trivial_ntt_matrix[i][j + row], in_consts[j + idx], Q);
          }
          else {
            plaintext_coefficients[j + col] += intel::hexl::MultiplyMod(trivial_ntt_matrix[i][j + row],
                                                                        omega * in_consts[j + idx], Q);
          }
        }
      }
    }
    rns_gsw_dot_product(acc[i], p, plaintext_coefficients, n, aut_ksk, priv_ksk);
  }

  for (size_t i = 0; i < n; i++) {
    gsw_half_RNSc_to_RNS_gsw(p[i], acc[i], priv_ksk);
  }
}

// Implementation of the part 2 over RGSW registers (that is using
// `rns_gsw_dot_product` for dot products.
// This function is slower than `intt_negacyclic_part2_quadratic_rnsc_rlwe`
// that compute the final part over RLWE registers, so it's only
// used for reference.
void intt_negacyclic_part2_quadratic(RNS_GSW * p, RNSc_half_GSW * acc,
                                     uint64_t * w, uint64_t n,
                                     uint64_t N_main,
                                     RNS_RLWE_KS_Key * aut_ksk, RNS_RLWE_KS_Key * priv_ksk) {
  const uint64_t Q = p[0]->samples[0]->b->Q;

  // Build the incomplete INTT matrix corresponding to the part 1 (quadratic) step of the INTT computation
  uint64_t trivial_ntt_matrix[N_main][N_main] = {0};
  for (size_t i = 0; i < N_main; i++) {
    uint64_t column_i[N_main];
    get_part2_intt_matrix_column(column_i, i, N_main, Q, INTT_QUADRATIC_LEVEL, NTT_INCOMPLETENESS);
    for (size_t j = 0; j < N_main; j++)
      trivial_ntt_matrix[j][i] = column_i[j];
  }

  for (size_t i = 0; i < n; i++) {
    uint64_t plaintext_coefficients[n] = {0};
    for (size_t j = 0; j < n; j++) {
      plaintext_coefficients[j] = trivial_ntt_matrix[i][j];
    }

    fprintf(stderr, "Part 2: %4ld / %4ld   ", i + 1, n);

    MEASURE_TIME("", 1, "rns_gsw_dot_product ", {
    perf_counters_start_if_measure_cache_misses_is_set();
    rns_gsw_dot_product(acc[i], p, plaintext_coefficients, n, aut_ksk, priv_ksk);
    perf_counters_stop_if_measure_cache_misses_is_set();
    })


    fflush(stderr);
  }
  fprintf(stderr, "\n");
  for (size_t i = 0; i < n; i++) {
    gsw_half_RNSc_to_RNS_gsw(p[i], acc[i], priv_ksk);
  }
}


// This function computes Part 2 of the INTT
void intt_negacyclic_part2_quadratic_rnsc_rlwe(RNSc_RLWE *out_tv, RNS_GSW * p,
                                     uint64_t * w, uint64_t n,
                                     uint64_t N_main,
                                     RNS_RLWE_KS_Key * aut_ksk, RNS_RLWE_KS_Key * priv_ksk) {
  const uint64_t Q = p[0]->samples[0]->b->Q;

  // Build the incomplete INTT matrix corresponding to the part 1 (quadratic) step of the INTT computation
  uint64_t trivial_ntt_matrix[N_main][N_main] = {0};

  for (size_t i = 0; i < N_main; i++) {
    uint64_t column_i[N_main];
    get_part2_intt_matrix_column(column_i, i, N_main, Q, INTT_QUADRATIC_LEVEL, NTT_INCOMPLETENESS);
    for (size_t j = 0; j < N_main; j++)
      trivial_ntt_matrix[j][i] = column_i[j];
  }

  #pragma omp parallel for
  for (size_t i = 0; i < n; i++) {
    uint64_t plaintext_coefficients[n] = {0};
    for (size_t j = 0; j < n; j++) {
      plaintext_coefficients[j] = trivial_ntt_matrix[i][j];
    }
    rnsc_rlwe_dot_product(out_tv[i], p, plaintext_coefficients, n, aut_ksk, priv_ksk);
  }
}


// This is the implementation of the fast INTT, that is not adequate for FHE because
// of the error growth. We do not use this implementation but kept it for future reference.
void _intt_rlwe_two_parts_iterative(RNSc_RLWE * tv, RNS_GSW * p, uint64_t * in_consts,
                                    RNSc_half_GSW * acc, uint64_t n,
                                    RNS_RLWE_KS_Key * aut_ksk,
                                    RNS_RLWE_KS_Key * priv_ksk){

  const uint64_t Q = p[0]->samples[0]->b->Q;

  uint64_t ws[n];
  uint64_t ws_inverse[n];
  uint64_t root_of_unity = intel::hexl::MinimalPrimitiveRoot((2 * n) >> NTT_INCOMPLETENESS, Q);
  calcRootsOfUnit(ws, root_of_unity, Q, n);
  calcRootsOfUnit(ws_inverse, intel::hexl::InverseMod(root_of_unity, Q), Q, n);

  fprintf(stderr, "Running trivial NTT (n = %ld, NTT_INCOMPLETENESS = %d).\n", n, NTT_INCOMPLETENESS);
  intt_negacyclic_part1_quadratic(p, in_consts, acc, ws, n, n, aut_ksk, priv_ksk);

  int k = ilog2(n);

  for (int i = INTT_QUADRATIC_LEVEL - 1; i >= 0; i--) {
    int d = 1 << (k - 1 - i);
    for (int j = 0; j < (1 << i); j++) {

      int phi_idx = bit_rev[bit_rev_index_for_n(n >> NTT_INCOMPLETENESS)][(1 << i) + j];
      uint64_t phi = ws_inverse[phi_idx];
      uint64_t inv2 = intel::hexl::InverseMod(2, Q);
      uint64_t coeffs0[2] = {inv2, inv2};
      uint64_t coeffs1[2] = {intel::hexl::MultiplyMod(inv2, phi, Q),
                             intel::hexl::MultiplyMod(inv2, intel::hexl::SubUIntMod(0, phi, Q), Q)};

      for (int u = (j * 2 * d); u < (j * 2 * d + d); u++) {
        RNS_GSW p_tmp[2] = {p[u], p[u + d]};

        rns_gsw_dot_product(acc[u], p_tmp, coeffs0, 2, aut_ksk, priv_ksk);
        rns_gsw_dot_product(acc[u + d], p_tmp, coeffs1, 2, aut_ksk, priv_ksk);

        gsw_half_RNSc_to_RNS_gsw(p[u], acc[u], priv_ksk);
        gsw_half_RNSc_to_RNS_gsw(p[u + d], acc[u + d], priv_ksk);

      }
    }
  }

  for (size_t k1 = 0; k1 < n; k1++){
    MACRO_rlwe_automorphism(tv[k1], tv[k1], 1);
    gsw_mul_RNSc_rlwe(tmp->samples[0], p[k1], tv[k1]);
    MACRO_rlwe_automorphism(tv[k1], tmp->samples[0], 1);
  }
}


void _intt_rlwe_two_parts_matrix(RNSc_RLWE * tv, RNS_GSW * p, uint64_t * in_consts,
                                 RNSc_half_GSW * acc, uint64_t n,
                                 RNS_RLWE_KS_Key * aut_ksk,
                                 RNS_RLWE_KS_Key * priv_ksk){

  const uint64_t Q = p[0]->samples[0]->b->Q;

  uint64_t ws[n];
  uint64_t ws_inverse[n];
  uint64_t root_of_unity = intel::hexl::MinimalPrimitiveRoot((2 * n) >> NTT_INCOMPLETENESS, Q);
  calcRootsOfUnit(ws, root_of_unity, Q, n);
  calcRootsOfUnit(ws_inverse, intel::hexl::InverseMod(root_of_unity, Q), Q, n);

  fprintf(stderr, "Part 1 of the INTT (n = %ld, NTT_INCOMPLETENESS = %d).\n", n, NTT_INCOMPLETENESS);
  fflush(stderr);
  MEASURE_TIME("", 1, "   -> ", {
  intt_negacyclic_part1_quadratic(p, in_consts, acc, ws, n, n, aut_ksk, priv_ksk);
  })

  fprintf(stderr, "Part 2 of the INTT\n");
  fflush(stderr);
  MEASURE_TIME("", 1, "   -> ", {
  intt_negacyclic_part2_quadratic_rnsc_rlwe(tv, p, ws, n, n, aut_ksk, priv_ksk);
  })
}


// Packed version of the blind rotate.
void packed_blind_rotate(RNSc_RLWE * tv, uint64_t * a_ntt, RNS_GSW * s, const uint64_t n, uint64_t root_of_unity_2nth, RNS_RLWE_KS_Key * aut_ksk, RNS_RLWE_KS_Key * priv_ksk){
  _intt_rlwe_two_parts_matrix(tv, s, a_ntt, acc, n, aut_ksk, priv_ksk);
}

RNS_GSW * new_rlwe_bootstrap_key(RNS_GSW_Key out_key, RLWE_Key in_key, uint64_t rou_2nth){
  assert(in_key->k == 1);
  const uint64_t N = in_key->N, Q = out_key->rlwe_key->Q, Q2 = in_key->q;
  uint64_t s_dft[N];
  array_additive_inverse_mod_switch(s_dft, in_key->s[0]->coeffs, Q2, Q, N);
  nc_ntt_forward(s_dft, s_dft, rou_2nth, Q, N);
  RNS_GSW * s = (RNS_GSW *) safe_malloc(sizeof(RNS_GSW)*N);
  for (size_t i = 0; i < N; i++){
    s[i] = (RNS_GSW) gsw_new_RNSc_sample(out_key, s_dft[i]);
    gsw_RNSc_to_RNS(s[i], (RNSc_GSW) s[i]);
  }
  return s;
} 

void packed_bootstrap_wo_extract(RNSc_RLWE * out, RNSc_RLWE * tv, RLWE in, uint64_t rou_2nth, RNS_GSW * s, RNS_RLWE_KS_Key * aut_ksk, RNS_RLWE_KS_Key * priv_ksk){
  assert(in->k == 1);
  const uint64_t n = in->a[0]->N, Q = out[0]->a->Q, Q2 = in->ntt->GetModulus(), l = s[0]->l;
  if(!acc) setup_tmp_buffer(n, Q, l, out[0]->a->ntt);
  array_mod_switch(in->b->coeffs, in->b->coeffs, Q2, Q, n);
  array_mod_switch(in->a[0]->coeffs, in->a[0]->coeffs, Q2, Q, n);
  for (size_t i = 0; i < n; i++){
    rlwe_RNSc_mul_by_xai(out[i], tv[i], (Q - in->b->coeffs[i]) % Q);
  }
  uint64_t a_ntt[n];
  nc_ntt_forward(a_ntt, in->a[0]->coeffs, rou_2nth, Q, n);
  packed_blind_rotate(out, a_ntt, s, n, rou_2nth, aut_ksk, priv_ksk);
}

// generates the complete keyset for the rlwe bootstrap and ks
RLWE_Bootstrap_KeySet new_rlwe_bootstrap_keyset(RLWE_Key in_key, RNS_GSW_Key gsw_key, LWE_Key lwe_key, uint64_t ks_l, uint64_t ks_base_bit){
  // parameters
  // ls and bases
  // const uint64_t ks_l = 23, ks_base_bit = 1;
  // paramters extracted from keys
  RNS_RLWE_Key rlwe_key = gsw_key->rlwe_key;
  const uint64_t Q = rlwe_key->Q, Ni = in_key->N, l = gsw_key->l;
  const uint64_t Q2 = in_key->q;
  intel::hexl::NTT * h_ntt = in_key->ntt;
  intel::hexl::NTT ** ntt = rlwe_key->s_RNS->ntt;
  const uint64_t p0 = ntt[0]->GetModulus();
  uint64_t rou_2nth = intel::hexl::MinimalPrimitiveRoot((2 * Ni) >> NTT_INCOMPLETENESS, Q);
  assert(rou_2nth != 0);
  // std::cout << "rou 2nth: " << rou_2nth << "\n";

  // generate RNS KS keys
  RNS_RLWE_KS_Key * priv_ksk = rlwe_new_RNS_priv_ks_key(rlwe_key, rlwe_key);
  RNS_RLWE_KS_Key * aut_ksk = rlwe_new_RNS_automorphism_keyset(rlwe_key);

  // generate RNS Bootstrapping key
  RNS_GSW * bk = new_rlwe_bootstrap_key(gsw_key, in_key, rou_2nth);

  // LWE ks key
  LWE_Key extracted_key = lwe_alloc_key(Q);
  rlwe_RNSc_extract_lwe_key(extracted_key, rlwe_key);
  extracted_key->q = Q2;
  LWE_KS_Key lwe_ksk = lwe_new_KS_key(lwe_key, extracted_key, ks_l, ks_base_bit);

  // packing key
  assert(in_key->q == lwe_key->q);
  RLWE_KS_Key packing_ksk = rlwe_new_full_packing_KS_key(in_key, lwe_key, ks_l, ks_base_bit);

  // create struct containing the keys
  RLWE_Bootstrap_KeySet res;
  res = (RLWE_Bootstrap_KeySet) safe_malloc(sizeof(*res));
  res->aut_ksk = aut_ksk;
  res->priv_ksk = priv_ksk;
  res->bk = bk;
  res->lwe_ksk = lwe_ksk;
  res->packing_ksk = packing_ksk;
  // allocate temporaries
  res->tmp_packed_in = rlwe_alloc_sample(Ni, 1, h_ntt);
  res->tmp_out = (RNSc_RLWE *) safe_malloc(sizeof(RNSc_RLWE)*Ni);
  for (size_t i = 0; i < Ni; i++) res->tmp_out[i] = (RNSc_RLWE) rlwe_new_RNS_trivial_sample_of_zero(Q, l, ntt);
  res->tmp_lwe_out = (LWE *) safe_malloc(sizeof(LWE)*Ni);
  for (size_t i = 0; i < Ni; i++) res->tmp_lwe_out[i] = lwe_alloc_sample(lwe_key->n, Q2);
  res->tmp_extracted_in = lwe_alloc_sample(Q, Q2); 
  res->s_tmp = (RNS_GSW *) safe_malloc(sizeof(RNS_GSW)*Ni);
  for (size_t i = 0; i < Ni; i++) res->s_tmp[i] = (RNS_GSW) gsw_new_RNSc_sample(gsw_key, 0);
  // save some parameters
  res->l = l;
  res->Q2 = Q2;
  res->Q = Q;
  res->rou_2nth = rou_2nth;
  return res;
}

void rlwe_bootstrap_and_ks(RLWE out, RLWE in, RNSc_RLWE tv, RLWE_Bootstrap_KeySet bks){
  const uint64_t Ni = in->b->N;
  RNSc_RLWE tv_vec[Ni];
  for (size_t i = 0; i < Ni; i++) tv_vec[i] = tv;

  // copy bk
  for (size_t i = 0; i < Ni; i++) gsw_RNS_copy(bks->s_tmp[i], bks->bk[i]);

  // run packed bootstrap
  packed_bootstrap_wo_extract(bks->tmp_out, tv_vec, in, bks->rou_2nth, bks->s_tmp, bks->aut_ksk, bks->priv_ksk);

  // mod switch to Q2, extract, and key switch to lwe_key
  for (size_t i = 0; i < Ni; i++){
    for (size_t j = 0; j < bks->l - 1 - SHRINK_LEVELS; j++){
      polynomial_base_reduce_RNSc_wo_free(bks->tmp_out[i]->a);
      polynomial_base_reduce_RNSc_wo_free(bks->tmp_out[i]->b);
    }
    rlwe_RNSc_mod_switch(bks->tmp_out[i], bks->Q2);
    rlwe_RNSc_extract_lwe(bks->tmp_extracted_in, bks->tmp_out[i], 0);
    lwe_keyswitch(bks->tmp_lwe_out[i], bks->tmp_extracted_in, bks->lwe_ksk);
    // reconstruct tmp_out buffer for the next boostrap
    bks->tmp_out[i]->a->l = bks->l;
    bks->tmp_out[i]->b->l = bks->l;
  }

  // packing
  rlwe_full_packing_keyswitch(out, bks->tmp_lwe_out, Ni, bks->packing_ksk);
}

void lwe_amortized_bootstrap(LWE * out, LWE * in, RNSc_RLWE tv, RLWE_Bootstrap_KeySet bks){
  const uint64_t Ni = bks->tmp_packed_in->b->N;

  fprintf(stderr, "[*] Initializing operation counters for lwe_amortized_bootstrap\n"); fflush(stderr);
  operation_counters_init();

  // packing
  rlwe_full_packing_keyswitch(bks->tmp_packed_in, in, Ni, bks->packing_ksk);

  // same TV for all lwes
  RNSc_RLWE tv_vec[Ni];
  for (size_t i = 0; i < Ni; i++) tv_vec[i] = tv;

  // copy bk
  for (size_t i = 0; i < Ni; i++) gsw_RNS_copy(bks->s_tmp[i], bks->bk[i]);

  // run packed bootstrap
  packed_bootstrap_wo_extract(bks->tmp_out, tv_vec, bks->tmp_packed_in, bks->rou_2nth, bks->s_tmp, bks->aut_ksk, bks->priv_ksk);

  // mod switch to Q2, extract, and key switch to lwe_key
  #pragma omp parallel for
  for (size_t i = 0; i < Ni; i++){
    for (size_t j = 0; j < bks->l - 1 - SHRINK_LEVELS; j++){
      polynomial_base_reduce_RNSc_wo_free(bks->tmp_out[i]->a);
      polynomial_base_reduce_RNSc_wo_free(bks->tmp_out[i]->b);
    }
    rlwe_RNSc_mod_switch(bks->tmp_out[i], bks->Q2);

    LWE tmp_extracted_in = lwe_alloc_sample(bks->Q, bks->Q2);
    rlwe_RNSc_extract_lwe(tmp_extracted_in, bks->tmp_out[i], 0);
    lwe_keyswitch(out[i], tmp_extracted_in, bks->lwe_ksk);
    free_lwe_sample(tmp_extracted_in);

    // reconstruct tmp_out buffer for the next boostrap
    bks->tmp_out[i]->a->l = bks->l;
    bks->tmp_out[i]->b->l = bks->l;
  }

  fprintf(stderr, "[*] Showing operation counters after lwe_amortized_bootstrap\n"); fflush(stderr);
  operation_counters_show();
}


void test_vector_from_array(RNSc_RLWE tv, uint64_t * in, uint64_t p, uint64_t size){
  const uint64_t p0 = tv->b->ntt[0]->GetModulus();
  uint64_t pi_prod = 1;
  for (size_t i = 1; i < tv->b->l; i++) pi_prod = intel::hexl::MultiplyMod(pi_prod, tv->b->ntt[i]->GetModulus(), p0);

  const double slot_size = (double) tv->b->Q/size;
  // printf("slot_size = %lf\n", slot_size);
  // const uint64_t slot_size = tv->b->Q/size;
  for (size_t i = 0; i < tv->b->Q; i++){
    uint64_t idx = round(i/slot_size);
    // const uint64_t in_ms = i/slot_size < size ? int_mod_switch(in[i/slot_size], p, p0) : 0;
    const uint64_t in_ms = i/slot_size < size ? int_mod_switch(in[idx], p, p0) : 0;
    tv->b->coeffs[0][i] = intel::hexl::MultiplyMod(in_ms, pi_prod, p0);
  }
}
