#include "rns-gsw.h"
#include "benchmark_util.h"
#include "omp.h"

// #define INTT_QUADRATIC_LEVEL 4
// #define NTT_INCOMPLETENESS 3

int INTT_QUADRATIC_LEVEL = 0;
int NTT_INCOMPLETENESS = 0;

int main(int argc, char const *argv[]) {
  // These are the original parameters:
  //    const uint64_t N = 1024, l = 3, bit_size = 49, p = 12289, p_start_size = 24, key_hw = 256;
  // ---

  if (argc == 2 && strcmp(argv[1], "print_header") == 0) {
    // Header
    printf("modulus_p,NTT_INCOMPLETENESS,INTT_QUADRATIC_LEVEL,N,l,bit_size,p_start_size,key_hw,N_prime,n_message_bits,std_err,dfr,OMP_NUM_THREADS,bootstrapping_time_ms\n");
    exit(0);
  }

  if (argc != 9) {
    printf("Usage: %s <prime_p> <n_message_bits> <N> <NTT_INCOMPLETENESS> <INTT_QUADRATIC_LEVEL> <l> <p_star_size> <key_hw>\n", argv[0]);
    exit(1);
  }


  uint64_t p = atoi(argv[1]);
  uint64_t n_message_bits = atoi(argv[2]);
  uint64_t N = atoi(argv[3]);
  NTT_INCOMPLETENESS = atoi(argv[4]);
  INTT_QUADRATIC_LEVEL = atoi(argv[5]);

  uint64_t l = atoi(argv[6]);
  uint64_t p_start_size = atoi(argv[7]);
  uint64_t key_hw = atoi(argv[8]);

  const uint64_t bit_size = 49;

  set_modulus_p(p);

  const uint64_t N_prime = next_power_of_2(p)<<1;

  // Generate primes
  auto p_star_list = intel::hexl::GeneratePrimes(1, p_start_size, true, N);
  const uint64_t p_star = p_star_list[0]; // p* is the output modulus
  auto Q = intel::hexl::GeneratePrimes(l, bit_size, true, N_prime);

  // Rings
  intel::hexl::NTT ** ntt = new_ntt_list(Q.data(), N_prime, l);
  intel::hexl::NTT * h_ntt = new intel::hexl::NTT(N, p_star);

  // generate keys
  // std::cout << "Generating private keys...\n";
  RLWE_Key in_key = rlwe_new_sparse_ternary_key(N, p_star, 1, key_hw, 1, h_ntt);
  RNS_RLWE_Key rlwe_key = rlwe_new_RNS_gaussian_key(p, l, 3.2, ntt, 1);
  RNS_GSW_Key gsw_key = gsw_new_RNS_key(rlwe_key);
  LWE_Key lwe_key = lwe_new_sparse_ternary_key(N, p_star, key_hw, 1);

  // RLWE_Key in_key = rlwe_new_sparse_ternary_key(N, p_star, 1, 0, 0, h_ntt);
  // RNS_RLWE_Key rlwe_key = rlwe_new_RNS_gaussian_key(p, l, 0, ntt, 0);
  // RNS_GSW_Key gsw_key = gsw_new_RNS_key(rlwe_key);
  // LWE_Key lwe_key = lwe_new_sparse_ternary_key(N, p_star, 0, 0);

  const uint64_t input_prec = n_message_bits,
                 input_base = 1<<input_prec;
                 // message_scale = p_star/input_base;
  const double message_scale = (double) p_star/input_base;


  // generate input
  LWE in[N];
  // Message is multiplied by the scaling.
  for (size_t i = 0; i < N; i++){
    in[i] = lwe_new_sample(round((i%input_base) * message_scale), lwe_key);
  }

  // generate test_vector
  uint64_t lut[input_base];
  for (size_t i = 0; i < input_base; i++){
    lut[i] = round(i*message_scale);
  }
  RNSc_RLWE tv = (RNSc_RLWE) rlwe_new_RNS_trivial_sample_of_zero(p, l, ntt);
  test_vector_from_array(tv, lut, p_star, input_base);


  // generate bootstrapping key
  fprintf(stderr, "Generate bootstrapping key... ");
  RLWE_Bootstrap_KeySet bks = new_rlwe_bootstrap_keyset(in_key, gsw_key, lwe_key, p_start_size, 1);
  fprintf(stderr, "Done.");
  fflush(stderr);

  // run rlwe bootstrap
  uint64_t start_bootstrap_time_usec = get_time();
  lwe_amortized_bootstrap(in, in, tv, bks);
  uint64_t finish_bootstrap_time_usec = get_time();
  double bootstrap_time_ms = (finish_bootstrap_time_usec - start_bootstrap_time_usec) / 1000.0;

  // decrypt
  uint64_t m_out[N];

  // round
  for (size_t i = 0; i < N; i++){
    uint64_t dec = round(lwe_phase(in[i], lwe_key) / message_scale);
    m_out[i] = dec % input_base;
  }
  print_array("Out:", m_out, N);

  double err = 0;
  for (size_t i = 0; i < N; i++){
    err += pow(mod_dist(m_out[i], i%input_base, input_base), 2);
  }

  double n_errors = 0;
  for (size_t i = 0; i < N; i++){
    n_errors += m_out[i] == (i % input_base) ? 0 : 1;
  }

  // printf("modulus_p,NTT_INCOMPLETENESS,INTT_QUADRATIC_LEVEL,N,l,bit_size,p_start_size,key_hw,N_prime,n_message_bits,std_err,dfr,OMP_NUM_THREADS,bootstrapping_time_ms\n");
  printf("%ld,", p);
  printf("%d,", NTT_INCOMPLETENESS);
  printf("%d,", INTT_QUADRATIC_LEVEL);
  printf("%ld,", N);
  printf("%ld,", l);
  printf("%ld,", bit_size);
  printf("%ld,", p_start_size);
  printf("%ld,", key_hw);
  printf("%ld,", N_prime);
  printf("%ld,", n_message_bits);
  printf("%lf,", sqrt(err / N));
  printf("%lf,", n_errors / N);
  printf("%d,", omp_get_max_threads());
  printf("%.4lf", bootstrap_time_ms);
  printf("\n");

  return 0;
}

