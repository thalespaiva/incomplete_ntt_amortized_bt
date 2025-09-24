#include "rns-gsw.h"

// #define DEBUG_NTT_PRINT // Prints inputs/outputs of NTT the computations

extern int NTT_INCOMPLETENESS;

void bit_rev_list(uint64_t * out, uint64_t * in, uint64_t n){
  // assert(n == 16 || n == 32 || n == 64 || n == 128 || n == 256 || n == 512 || n == 1024 || n == 2048);
  for (size_t i = 0; i < n; i++){
    out[i] = in[bit_rev[bit_rev_index_for_n(n)][i]];
  }
}

void ntt_forward(uint64_t * out, uint64_t * in, uint64_t * ws, uint64_t Q, uint64_t n){
  assert(in != out);
  bit_rev_list(out, in, n);
  uint64_t m = 2;
  while(m <= n){
    const uint64_t w_m = ws[n/m];
    uint64_t w = 1;
    for (size_t j = 0; j < m/2; j++){
      // const uint64_t w = ws[j*n/m];
      for (size_t k = 0; k < n-1; k+=m){
        const uint64_t t = intel::hexl::MultiplyMod(w, out[k + j + m/2], Q);
        const uint64_t u = out[k + j];
        out[k + j] = intel::hexl::AddUIntMod(u, t, Q);
        out[k + j + m/2] = intel::hexl::SubUIntMod(u, t, Q);
      }
      w = intel::hexl::MultiplyMod(w, w_m, Q);
    }
    m *= 2;
  }
}

void calcRootsOfUnit(uint64_t * out, uint64_t min_RoU, uint64_t q, uint64_t n){
  uint64_t r = 1;
  for (size_t i = 0; i < n; i++){
    out[i] = r;
    r = intel::hexl::MultiplyMod(r, min_RoU, q);
  }
}

void entry_wise_prod(uint64_t * out, uint64_t * in1, uint64_t * in2, uint64_t n){
  for (size_t i = 0; i < n; i++){
    out[i] = in1[i] * in2[i];
  }
}

// negacyclic forward ntt
void nc_ntt_forward_old(uint64_t * out, uint64_t * in, uint64_t rou_2nth, uint64_t q, uint64_t n){
  uint64_t rou_vec[n], in2[n];
  calcRootsOfUnit(rou_vec, rou_2nth, q, n);
  entry_wise_prod(in2, in, rou_vec, n);
  const uint64_t rou = intel::hexl::MultiplyMod(rou_2nth, rou_2nth, q);
  calcRootsOfUnit(rou_vec, rou, q, n);
  ntt_forward(out, in2, rou_vec, q, n);
}


void nc_ntt_forward(uint64_t * out, uint64_t * in, uint64_t rou_2nth,
                   uint64_t Q, uint64_t n) {
  uint64_t ws[n];
  calcRootsOfUnit(ws, rou_2nth, Q, n);
  uint64_t ntt_in[n];
  for (size_t i = 0; i < n; i++) {
    ntt_in[i] = in[i];
  }

#ifdef DEBUG_NTT_PRINT
  fprintf(stderr, "[nc_ntt_forward] in: ");
  for (size_t i = 0; i < n; i++) {
    fprintf(stderr, "%ld, ", ntt_in[i]);
  }
  fprintf(stderr, "\n");
#endif

  int n_levels = ilog2(n) - NTT_INCOMPLETENESS;
  int k = ilog2(n);
  for (int i = 0; i < n_levels; i++) {
    int d = 1 << (k - 1 - i);
    for (int j = 0; j < (1 << i); j++) {
      int phi_idx = bit_rev[bit_rev_index_for_n(n >> NTT_INCOMPLETENESS)][(1 << i) + j];
      uint64_t phi = ws[phi_idx];
      for (int u = (j * 2 * d); u < (j * 2 * d + d); u++) {
        uint64_t t0 = ntt_in[u];
        uint64_t t1 = intel::hexl::MultiplyMod(phi, ntt_in[u + d], Q);
        ntt_in[u] = intel::hexl::AddUIntMod(t0, t1, Q);
        ntt_in[u + d] = intel::hexl::SubUIntMod(t0, t1, Q);
      }
    }
  }
  // Uncomment the line below to get the entries in bit-reversed order
  // bit_rev_list(out, ntt_in, n);
  for (size_t i = 0; i < n; i++) {
    out[i] = ntt_in[i];
  }

#ifdef DEBUG_NTT_PRINT
  fprintf(stderr, "out: ");
  for (size_t i = 0; i < n; i++) {
    fprintf(stderr, "%ld, ", out[i]);
  }
  fprintf(stderr, "\n");
#endif

}


void nc_ntt_forward_leveled(uint64_t * out, uint64_t * in, uint64_t rou, size_t rou_order,
                            uint64_t Q, uint64_t n, int incompleteness) {
  uint64_t ws[n];
  calcRootsOfUnit(ws, rou, Q, n);
  uint64_t ntt_in[n];
  for (size_t i = 0; i < n; i++) {
    ntt_in[i] = in[i];
  }

#ifdef DEBUG_NTT_PRINT
  fprintf(stderr, "in: ");
  for (size_t i = 0; i < n; i++) {
    fprintf(stderr, "%ld, ", ntt_in[i]);
  }
  fprintf(stderr, "\n");
#endif

  int bit_rev_length = rou_order/2;

  int n_levels = ilog2(n) - incompleteness;
  int k = ilog2(n);
  for (int i = 0; i < n_levels; i++) {
    int d = 1 << (k - 1 - i);
    for (int j = 0; j < (1 << i); j++) {
      int phi_idx = bit_rev[bit_rev_index_for_n(bit_rev_length)][(1 << i) + j];
      uint64_t phi = ws[phi_idx];
      for (int u = (j * 2 * d); u < (j * 2 * d + d); u++) {
        uint64_t t0 = ntt_in[u];
        uint64_t t1 = intel::hexl::MultiplyMod(phi, ntt_in[u + d], Q);
        ntt_in[u] = intel::hexl::AddUIntMod(t0, t1, Q);
        ntt_in[u + d] = intel::hexl::SubUIntMod(t0, t1, Q);
      }
    }
  }
  // Uncomment the line below to get the entries in bit-reversed order
  // bit_rev_list(out, ntt_in, n);
  for (size_t i = 0; i < n; i++) {
    out[i] = ntt_in[i];
  }

#ifdef DEBUG_NTT_PRINT
  fprintf(stderr, "out: ");
  for (size_t i = 0; i < n; i++) {
    fprintf(stderr, "%ld, ", out[i]);
  }
  fprintf(stderr, "\n");
#endif

}


// Computes the inverse incomplete NTT with a parameterized incompleteness level.
void nc_ntt_inverse_leveled(uint64_t * out, uint64_t * in, uint64_t rou, size_t rou_order,
                            uint64_t Q, uint64_t n, int incompleteness) {
  uint64_t ws[n];
  calcRootsOfUnit(ws, rou, Q, n);
  uint64_t ntt_in[n];
  for (size_t i = 0; i < n; i++) {
    ntt_in[i] = in[i];
  }

#ifdef DEBUG_NTT_PRINT
  fprintf(stderr, "[nc_ntt_inverse_leveled] in: ");
  for (size_t i = 0; i < n; i++) {
    fprintf(stderr, "%ld, ", ntt_in[i]);
  }
  fprintf(stderr, "\n");
#endif

  int bit_rev_length = rou_order/2;
  int n_levels = ilog2(n) - incompleteness;
  int k = ilog2(n);

  for (int i = n_levels - 1; i >= 0; i--) {
    int d = 1 << (k - 1 - i);
    for (int j = 0; j < (1 << i); j++) {

      int phi_idx = bit_rev[bit_rev_index_for_n(bit_rev_length)][(1 << i) + j];

      uint64_t phi = intel::hexl::InverseMod(ws[phi_idx], Q);

      for (int u = (j * 2 * d); u < (j * 2 * d + d); u++) {
        uint64_t t0 = intel::hexl::AddUIntMod(ntt_in[u], ntt_in[u + d], Q);
        uint64_t t1 = intel::hexl::SubUIntMod(ntt_in[u], ntt_in[u + d], Q);
        ntt_in[u] = t0;
        ntt_in[u + d] = intel::hexl::MultiplyMod(phi, t1, Q);
      }
    }
  }

  // Uncomment the line below to get the entries in bit-reversed order
  // bit_rev_list(out, ntt_in, n);
  uint64_t norm_factor = intel::hexl::InverseMod(n >> incompleteness, Q);
  for (size_t i = 0; i < n; i++) {
    out[i] = intel::hexl::MultiplyMod(ntt_in[i], norm_factor, Q);
  }

#ifdef DEBUG_NTT_PRINT
  fprintf(stderr, "[nc_ntt_inverse_leveled] out: ");
  for (size_t i = 0; i < n; i++) {
    fprintf(stderr, "%ld, ", out[i]);
  }
  fprintf(stderr, "\n");
#endif

}
