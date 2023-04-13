#include "packing.h"
#include "params.h"
#include "polydbl.h"
#include "randombytes.h"
#include "sampler.h"
#include "sign.h"
#include "symmetric.h"
#include <stdio.h>

#define REP 10000

int main() {
    uint8_t seed[SEEDBYTES];
    uint8_t buf[STREAM128_BLOCKBYTES];
    uint16_t nonce = 0;
    randombytes(seed, sizeof seed);
    FILE *f;
    f = fopen("gaussian.txt", "w");
    fprintf(f, "%d\n", REP);

    for (int i = 0; i < REP; ++i) {
        double temp[1];
        stream128_state state;
        stream128_init(&state, seed, nonce++);
        stream128_squeezeblocks(buf, 1, &state);
        if (sampler_gaussian(temp, 1, buf, STREAM128_BLOCKBYTES) == 0) {
            printf("not sampled!!!\n");
        }
        fprintf(f, "%.16e\n", temp[0]);
    }

    fclose(f);
    f = fopen("uniformdbl.txt", "w");
    fprintf(f, "%d\n", REP);

    for (int i = 0; i < REP; ++i) {
        stream128_state state;
        stream128_init(&state, seed, nonce++);
        stream128_squeezeblocks(buf, 1, &state);
        fprintf(f, "%.16e\n", sampler_uniformdbl(buf));
    }

    fclose(f);
    f = fopen("hyperball_radius.txt", "w");
    fprintf(f, "%d\n", REP);

    for (int i = 0; i < REP; ++i) {
        polydblvecl z0;
        polydblvecl z1;
        polydblveclk_uniform_hyperball(&z0, &z1, seed, nonce++);
        fprintf(f, "%.16e\n", polydblveclk_norm2(&z0, &z1));
    }

    fclose(f);
    // f = fopen("z_radius.txt", "w");
    // fprintf(f, "%d\n", REP);

    // {
    //     uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    //     uint8_t sk[CRYPTO_SECRETKEYBYTES];
    //     uint8_t sig[CRYPTO_BYTES];
    //     size_t siglen;
    //     polyvecm z;
    //     uint8_t c[SEEDBYTES];
    //     crypto_sign_keypair(pk, sk);

    //     for (int i = 0; i < REP; ++i) {
    //         crypto_sign_signature(sig, &siglen, pk, CRHBYTES, sk);
    //         unpack_sig(c, &z, sig);
    //         fprintf(f, "%.16e\n", polyvecm_norm2(&z));
    //     }
    // }

    // fclose(f);

    return 0;
}