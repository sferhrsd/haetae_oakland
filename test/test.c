#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "poly.h"
#include "ntt.h"
#include "randombytes.h"

#define NTESTS 100

// NTT^-1(NTT(x) = x
void testNTT() {
    uint8_t seed[SEEDBYTES];
    uint16_t nonce = 0;
    randombytes(seed, SEEDBYTES);

    for (int t = 0; t < NTESTS; t++) {
        poly a;
        poly b;
        poly_uniform(&a, seed, nonce++);
        b = a;

        ntt(&a.coeffs[0]);
        invntt_tomont(&a.coeffs[0]);     
        for (int i = 0; i < N; ++i) {
            a.coeffs[i] = montgomery_reduce((int64_t)a.coeffs[i]);
        }

        poly_freeze(&a);
        poly_freeze(&b);

        int equal = !memcmp(&a, &b, sizeof(poly));
        assert(equal && "NTT^-1(NTT(x) != x");
    }
}

void testNTTmul() {
    uint8_t seed[SEEDBYTES];
    uint16_t nonce = 0;
    randombytes(seed, SEEDBYTES);

    for (int t = 0; t < NTESTS; t++) {
        poly a;
        poly b;
        poly c;
        poly_uniform(&a, seed, nonce++);
        poly_uniform(&b, seed, nonce++);

        poly_naivemul(&c, &a, &b);

        ntt(&a.coeffs[0]);
        ntt(&b.coeffs[0]);
        poly_pointwise_montgomery(&a, &a, &b);
        invntt_tomont(&a.coeffs[0]);     

        poly_freeze(&a);
        poly_freeze(&c);

        int equal = !memcmp(&a, &c, sizeof(poly));
        assert(equal && "NTT Mul != Schoolbock Mul");
    }
}


int main(void) {
    testNTT();
    testNTTmul();
    printf("TESTS PASSED\n");
    return 0;
}