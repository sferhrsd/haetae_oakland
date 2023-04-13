#ifndef SYMMETRIC_H
#define SYMMETRIC_H

#include "params.h"
#include <stdint.h>
#include "fips202.h"

// Cryptographic XOF function: blake3 or shake256
#ifdef HAETAE_USE_BLAKE // XOF: use blake

#include "blake3.h"

typedef blake3_hasher xof256_state;

#define haetae_blake_absorbe_once                                               \
    HAETAE_NAMESPACE(haetae_blake_absorbe_once)
void haetae_blake_absorbe_once(blake3_hasher *state, const uint8_t *in, size_t inlen);

#define haetae_blake_absorbe_twice                                              \
    HAETAE_NAMESPACE(haetae_blake_absorbe_twice)
void haetae_blake_absorbe_twice(blake3_hasher *state, const uint8_t *in1, size_t in1len, 
                                const uint8_t *in2, size_t in2len);

#define XOF256_BLOCKBYTES BLAKE3_OUT_LEN

#define xof256_absorbe_once(STATE, IN, IN_LEN)                                  \
    haetae_blake_absorbe_once(STATE, IN, IN_LEN)
#define xof256_absorbe_twice(STATE, IN1, IN1_LEN, IN2, IN2_LEN)                 \
    haetae_blake_absorbe_twice(STATE, IN1, IN1_LEN, IN2, IN2_LEN)
#define xof256_squeeze(OUT, OUT_LEN, STATE)                                     \
    blake3_hasher_finalize(STATE, OUT, OUT_LEN) 
#define xof256_squeezeblocks(OUT, OUTBLOCKS, STATE)                             \
    blake3_hasher_finalize(STATE, OUT, OUTBLOCKS*BLAKE3_OUT_LEN)

#else // XOF: use shake256

typedef keccak_state xof256_state;

#define haetae_shake256_absorb_twice                                           \
    HAETAE_NAMESPACE(haetae_shake256_absorb_twice)
void haetae_shake256_absorb_twice(keccak_state *state, const uint8_t *in1,
                                size_t in1len, const uint8_t *in2, size_t in2len);

#define XOF256_BLOCKBYTES SHAKE256_RATE

#define xof256_absorbe_once(STATE, IN, IN_LEN)                                  \
    shake256_absorb_once(STATE, IN, IN_LEN)
#define xof256_absorbe_twice(STATE, IN, IN_LEN, IN2, IN2_LEN)                   \
    haetae_shake256_absorb_twice(STATE, IN, IN_LEN, IN2, IN2_LEN) 
#define xof256_squeeze(OUT, OUT_LEN, STATE)                                     \
    shake256_squeeze(OUT, OUT_LEN, STATE)
#define xof256_squeezeblocks(OUT, OUTBLOCKS, STATE)                             \
    shake256_squeezeblocks(OUT, OUTBLOCKS, STATE)

#endif // XOF


// Stream function: aes256 or blake3 or shake128|256
#ifdef HAETAE_USE_AES // stream: aes256

#include "aes256ctr.h"

typedef aes256ctr_ctx stream128_state;
typedef aes256ctr_ctx stream256_state;

#define HAETAE_aes256ctr_init HAETAE_NAMESPACE(HAETAE_aes256ctr_init)
void HAETAE_aes256ctr_init(aes256ctr_ctx *state,
                              const uint8_t key[32],
                              uint16_t nonce);

#define STREAM128_BLOCKBYTES AES256CTR_BLOCKBYTES
#define STREAM256_BLOCKBYTES AES256CTR_BLOCKBYTES

#define stream128_init(STATE, SEED, NONCE) \
        HAETAE_aes256ctr_init(STATE, SEED, NONCE)
#define stream128_squeezeblocks(OUT, OUTBLOCKS, STATE) \
        aes256ctr_squeezeblocks(OUT, OUTBLOCKS, STATE)
#define stream256_init(STATE, SEED, NONCE) \
        HAETAE_aes256ctr_init(STATE, SEED, NONCE)
#define stream256_squeezeblocks(OUT, OUTBLOCKS, STATE) \
        aes256ctr_squeezeblocks(OUT, OUTBLOCKS, STATE)

#elif defined (HAETAE_USE_BLAKE)   // stream: blake3

typedef blake3_hasher stream128_state;
typedef blake3_hasher stream256_state;

#define haetae_blake_stream_init                                                \
    HAETAE_NAMESPACE(haetae_blake_stream_init)
void haetae_blake128_stream_init(blake3_hasher *state,
                                 const uint8_t seed[SEEDBYTES], uint16_t nonce);

#define haetae_blake_stream_init                                                \
    HAETAE_NAMESPACE(haetae_blake_stream_init)
void haetae_blake256_stream_init(blake3_hasher *state,
                                 const uint8_t seed[CRHBYTES], uint16_t nonce);

#define STREAM128_BLOCKBYTES BLAKE3_OUT_LEN
#define STREAM256_BLOCKBYTES BLAKE3_OUT_LEN

// for blake the is actually no differentiation in security level 
// but we need two differen versions of stream_init anyway:
#define stream128_init(STATE, SEED, NONCE)                                     \
    haetae_blake128_stream_init(STATE, SEED, NONCE)
#define stream128_squeezeblocks(OUT, OUTBLOCKS, STATE)                         \
    blake3_hasher_finalize(STATE, OUT, OUTBLOCKS*BLAKE3_OUT_LEN)
#define stream256_init(STATE, SEED, NONCE)                                     \
    haetae_blake256_stream_init(STATE, SEED, NONCE)
#define stream256_squeezeblocks(OUT, OUTBLOCKS, STATE)                         \
    blake3_hasher_finalize(STATE, OUT, OUTBLOCKS*BLAKE3_OUT_LEN)

#else // stream: shake128 and shake256

typedef keccak_state stream128_state;
typedef keccak_state stream256_state;

#define haetae_shake128_stream_init                                            \
    HAETAE_NAMESPACE(haetae_shake128_stream_init)
void haetae_shake128_stream_init(keccak_state *state,
                                 const uint8_t seed[SEEDBYTES], uint16_t nonce);

#define haetae_shake256_stream_init                                            \
    HAETAE_NAMESPACE(haetae_shake256_stream_init)
void haetae_shake256_stream_init(keccak_state *state,
                                 const uint8_t seed[CRHBYTES], uint16_t nonce);

#define STREAM128_BLOCKBYTES SHAKE128_RATE
#define STREAM256_BLOCKBYTES SHAKE256_RATE

#define stream128_init(STATE, SEED, NONCE)                                     \
    haetae_shake128_stream_init(STATE, SEED, NONCE)
#define stream128_squeezeblocks(OUT, OUTBLOCKS, STATE)                         \
    shake128_squeezeblocks(OUT, OUTBLOCKS, STATE)
#define stream256_init(STATE, SEED, NONCE)                                     \
    haetae_shake256_stream_init(STATE, SEED, NONCE)
#define stream256_squeezeblocks(OUT, OUTBLOCKS, STATE)                         \
    shake256_squeezeblocks(OUT, OUTBLOCKS, STATE)

#endif // stream

#endif //SYMMETRIC_H
