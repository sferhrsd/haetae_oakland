// clang-format off
#ifndef HAETAE_PARAMS_H
#define HAETAE_PARAMS_H

#include "config.h"

#define SEEDBYTES 32
#define CRHBYTES 64
#define N 256
#define ROOT_OF_UNITY 3

#define Q 64513
#define DQ (Q << 1)// 2Q

#define BETA 256
#define HALF_BETA (BETA >> 1)
#define LOG_BETA 8

#define MAX_ENC_H_BYTES (2*N*K*2+4) // TODO CHECK
#define MAX_ENC_HB_Z1_BYTES (N*L*2+4) // size_of(Encoded(HighBits(z1))) <= N * L * SCALE_BITS / 8bits + 4

#if HAETAE_MODE == 2
#define K 2
#define L 4
#define ETA 1
#define TAU 58
#define B0 9388.96
#define B1 9382.26
#define B2 12320.79
#define GAMMA 48.83
#define LN 8192 //6161 // Large N
#define SQNM 39.191835884530846 // \sqrt(n * m)
#define D 1
#define CRYPTO_BYTES 1463

#define ALPHA 512
#define LOG_ALPHA 9

#define POLYB1_PACKEDBYTES 480     // 15bits * N / 8bits

#elif HAETAE_MODE == 3
#define K 3
#define L 6
#define ETA 1
#define TAU 80
#define B0 17773.20
#define B1 17766.15
#define B2 21365.10
#define GAMMA 57.68
#define LN 8192//7045 // Large N
#define SQNM 48.0
#define D 1
#define CRYPTO_BYTES 2337

#define ALPHA 512
#define LOG_ALPHA 9

#define POLYB1_PACKEDBYTES 480     // 15bits * N / 8bits

#elif HAETAE_MODE == 5
#define K 4
#define L 7
#define ETA 1
#define TAU 60
#define B0 20614.9815
#define B1 20609.9152
#define B2 23740.4482
#define GAMMA 55.13
#define LN 8192//7254 // Large N
#define SQNM 53.0659966456864
#define D 0
#define CRYPTO_BYTES 2908

#define ALPHA 256
#define LOG_ALPHA 8

#define POLYB1_PACKEDBYTES 512     // 16bits * N / 8bits

#endif



#define HALF_ALPHA (ALPHA >> 1) // ALPHA / 2

#define B0SQ ((uint64_t)(B0*B0))
#define B1SQ ((uint64_t)(B1*B1))
#define B2SQ ((uint64_t)(B2*B2))

#define POLYQ_PACKEDBYTES 512  // 16bits * N / 8bits
#define POLY2Q_PACKEDBYTES 544 // 17bits * N / 8bits


#define M (L-1)

#if ETA == 1
#define POLYETA_PACKEDBYTES 64
#define POLY2ETA_PACKEDBYTES 96 
#elif ETA == 2
#define POLYETA_PACKEDBYTES 96
#endif

#define POLYC_PACKEDBYTES 32       // 1bit * N / 8bits
#define POLY_HIGHBITS_PACKEDBYTES (N * 9 / 8)
#define POLYVECK_HIGHBITS_PACKEDBYTES (POLY_HIGHBITS_PACKEDBYTES * K)
#define POLYVECK_BYTES (K * N * sizeof(int32_t))
#define POLYVECL_BYTES (L * N * sizeof(int32_t))

#define CRYPTO_PUBLICKEYBYTES (SEEDBYTES + K * POLYQ_PACKEDBYTES)                                      // seed + b
#if D == 1
#define CRYPTO_SECRETKEYBYTES (CRYPTO_PUBLICKEYBYTES + M * POLYETA_PACKEDBYTES + K * POLY2ETA_PACKEDBYTES + SEEDBYTES)  // pk + s + K
#elif D == 0
#define CRYPTO_SECRETKEYBYTES (CRYPTO_PUBLICKEYBYTES + (M + K) * POLYETA_PACKEDBYTES + SEEDBYTES)  // pk + s + K
#else
#error
#endif
#endif
// clang-format on
