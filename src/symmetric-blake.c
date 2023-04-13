#include "blake3.h"
#include "params.h"
#include "symmetric.h"
#include <stdint.h>

void haetae_blake128_stream_init(blake3_hasher *state,
                                 const uint8_t seed[SEEDBYTES], uint16_t nonce) {
    uint8_t t[2];
    t[0] = nonce;
    t[1] = nonce >> 8;

    blake3_hasher_init(state);
    blake3_hasher_update(state, seed, SEEDBYTES);
    blake3_hasher_update(state, t, 2);
}

void haetae_blake256_stream_init(blake3_hasher *state,
                                 const uint8_t seed[CRHBYTES], uint16_t nonce) {
    uint8_t t[2];
    t[0] = nonce;
    t[1] = nonce >> 8;

    blake3_hasher_init(state);
    blake3_hasher_update(state, seed, CRHBYTES);
    blake3_hasher_update(state, t, 2);
}

void haetae_blake_absorbe_once(blake3_hasher *state, const uint8_t *in, size_t inlen) {
    blake3_hasher_init(state);
    blake3_hasher_update(state, in, inlen);
}

void haetae_blake_absorbe_twice(blake3_hasher *state, const uint8_t *in1, size_t in1len, 
                                const uint8_t *in2, size_t in2len) {
    blake3_hasher_init(state);
    blake3_hasher_update(state, in1, in1len);
    blake3_hasher_update(state, in2, in2len);
}
