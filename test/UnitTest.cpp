/*
  RLizard
  Copyright (C) 2021-present CryptoLab, Inc.
  This software is distributed under the dual license.
  The two licenses are:
  1. GNU Affero General Public License v3.0 (LICENSE-AGPL.txt)
  2. CryptoLab's Commercial License
  You can choose either of the two licenses at your option.
  To lift the constraints from GNU Affero General Public License, please
  contact us via cryptolab@cryptolab.co.kr.
  Licenses for Third-Party Components can be found in LICENSE-3rd-party.txt.
 */

#include "gtest/gtest.h"
#include <string>

extern "C" {
#include "ntt.h"
#include "packing.h"
#include "params.h"
#include "poly.h"
#include "polymat.h"
#include "polyvec.h"
#include "randombytes.h"
#include "sign.h"
}

#define NTESTS 1000
#define PRINT_FLAG 0

template <class Ty>
static bool arrayEq(const Ty *a, const Ty *b, unsigned len) {
    for (unsigned i = 0; i < len; ++i) {
        if (a[i] != b[i])
            return false;
    }
    return true;
}

static bool polyEqDQ(const poly *a, const poly *b) {
    for (unsigned j = 0; j < N; ++j) {
        if ((a->coeffs[j] - b->coeffs[j]) % DQ)
            return false;
    }
    return true;
}

static bool polyEqQ(const poly *a, const poly *b) {
    for (unsigned j = 0; j < N; ++j) {
        if ((a->coeffs[j] - b->coeffs[j]) % Q)
            return false;
    }
    return true;
}

static bool isSigDiff(const uint8_t sig1[CRYPTO_BYTES],
                      const uint8_t sig2[CRYPTO_BYTES]) {

    unsigned int i;
    polyvecm z1, z2;
    uint8_t c1[SEEDBYTES], c2[SEEDBYTES];

    unpack_sig(c1, &z1, sig1);
    unpack_sig(c2, &z2, sig2);
    for (i = 0; i < M; ++i) {
        if (!polyEqDQ(&z1.vec[i], &z2.vec[i]))
            return true;
    }
    if (memcmp(c1, c2, SEEDBYTES))
        return true;
    return false;
}

#define MUTATE(arr, len)                                                       \
    if (rand() % 2 == 0) {                                                     \
        int picked = rand();                                                   \
        if (picked != 0) {                                                     \
            (arr)[rand() % (len)] += picked;                                   \
        }                                                                      \
    }

void mutateSignature(uint8_t sig[CRYPTO_BYTES]) {
    uint8_t sig_org[CRYPTO_BYTES];
    memcpy(sig_org, sig, CRYPTO_BYTES);
    while (!isSigDiff(sig, sig_org)) {
        MUTATE(sig, CRYPTO_BYTES);
    }
}
#undef MUTATE

static void polyq_naivemul(poly *c, const poly *a, const poly *b) {
    unsigned int i, j;
    int32_t r[2 * N] = {0};

    for (i = 0; i < N; ++i)
        for (j = 0; j < N; ++j)
            r[i + j] = (r[i + j] + (int64_t)a->coeffs[i] * b->coeffs[j]) % Q;

    for (i = N; i < 2 * N; ++i)
        r[i - N] = (r[i - N] - r[i]) % Q;

    for (i = 0; i < N; ++i)
        c->coeffs[i] = r[i];
}

void testPolyNTT() {
    uint8_t seed[SEEDBYTES];
    uint16_t nonce = 0;
    randombytes(seed, sizeof(seed));
    for (int t = 0; t < NTESTS; t++) {
        poly a;
        poly b;
        poly_uniform(&a, seed, nonce++);
        for (int i = 0; i < N; ++i)
            b.coeffs[i] = a.coeffs[i];
        poly_ntt(&a);
        for (int i = 0; i < N; ++i) {
            a.coeffs[i] = (int64_t)a.coeffs[i] * (-432) % Q;
        }
        poly_invntt_tomont(&a);
        ASSERT_TRUE(polyEqQ(&a, &b));
    }
}

void testPolyMult() {
    uint8_t seed[SEEDBYTES];
    uint16_t nonce = 0;
    randombytes(seed, sizeof(seed));
    for (int t = 0; t < NTESTS; t++) {
        poly a;
        poly b;
        poly c;
        poly d;
        poly_uniform(&a, seed, nonce++);
        poly_uniform(&b, seed, nonce++);
        polyq_naivemul(&d, &a, &b);
        poly_ntt(&a);
        poly_ntt(&b);
        poly_pointwise_montgomery(&c, &a, &b);
        poly_invntt_tomont(&c);
        ASSERT_TRUE(polyEqQ(&c, &d));
    }
}

void testPolyInv() {
    uint8_t seed[SEEDBYTES];
    uint16_t nonce = 0;
    randombytes(seed, sizeof(seed));
    for (int t = 0; t < NTESTS; t++) {
        poly a;
        poly b;
        poly c;
        do {
            poly_uniform(&a, seed, nonce++);
            poly_ntt(&a);
        } while (!poly_pointwise_inverse_tomontsq(&b, &a));

        // checks if a * a^(-1) == MONTSQ(before ntt)
        for (int i = 0; i < N; ++i) {
            ASSERT_EQ((a.coeffs[i] * b.coeffs[i] - MONTSQ) % Q, 0);
        }
        poly_pointwise_montgomery(&c, &a, &b);
        poly_invntt_tomont(&c);
        for (int i = 0; i < N; ++i) { // checks if a * a^(-1) == MONTSQ
            int32_t val = i == 0 ? MONTSQ : 0;
            ASSERT_EQ((c.coeffs[i] - val) % Q, 0);
        }
    }
}

void testPolyGen() {
    uint8_t seed[CRHBYTES];
    uint16_t nonce = 0;
    for (int t = 0; t < NTESTS; t++) {
        randombytes(seed, sizeof(seed));
        poly a;
        poly_uniform(&a, seed, nonce++);
        if (PRINT_FLAG) {
            for (int i = 0; i < N; ++i) {
                printf("%d ", a.coeffs[i]);
            }
            printf("\n");
        }
        for (int i = 0; i < N; ++i) {
            ASSERT_TRUE(a.coeffs[i] >= 0 && a.coeffs[i] < Q);
        }
        poly_uniform_eta(&a, seed, nonce++);
        if (PRINT_FLAG) {
            for (int i = 0; i < N; ++i) {
                printf("%d ", a.coeffs[i]);
            }
            printf("\n");
        }
        for (int i = 0; i < N; ++i) {
            ASSERT_TRUE(a.coeffs[i] >= -ETA && a.coeffs[i] <= ETA);
        }
        poly_challenge(&a, seed);
        int onecnt = 0;
        for (int i = 0; i < N; ++i) {
            ASSERT_TRUE(0 <= a.coeffs[i] && a.coeffs[i] <= 1);
            onecnt += a.coeffs[i];
        }
        ASSERT_EQ(onecnt, TAU);
    }
}

void testPolyFreeze() {
    for (int t = 0; t < NTESTS; t++) {
        poly a;
        poly b;
        for (int i = 0; i < N; ++i) {
            b.coeffs[i] = a.coeffs[i] = rand();
        }
        poly_freeze2q(&a);
        for (int i = 0; i < N; ++i) {
            ASSERT_TRUE(a.coeffs[i] >= 0 && a.coeffs[i] <= DQ);
            ASSERT_EQ((a.coeffs[i] - b.coeffs[i]) % Q, 0);
        }
    }
}

void testPolyCRT() {
    uint8_t seed[SEEDBYTES];
    uint16_t nonce = 0;
    randombytes(seed, sizeof(seed));
    for (int t = 0; t < NTESTS; t++) {
        poly a, b, c;
        poly_uniform(&a, seed, nonce++);
        poly_uniform(&b, seed, nonce++);
        for (int i = 0; i < N; ++i) {
            b.coeffs[i] = b.coeffs[i] % 2; // sending to mod 2 space
        }
        c = a;
        poly_crt(&c, &b);
        for (int i = 0; i < N; ++i) {
            ASSERT_EQ((c.coeffs[i] - a.coeffs[i]) % Q, 0);
            ASSERT_EQ((c.coeffs[i] - b.coeffs[i]) % 2, 0);
        }
    }
}

void testPolyPack() {
    uint8_t seed[CRHBYTES];
    uint16_t nonce = 0;
    randombytes(seed, sizeof(seed));
    poly a;
    poly b;
    uint8_t rq[POLYQ_PACKEDBYTES];
    uint8_t reta[POLYETA_PACKEDBYTES];
    uint8_t r2q[POLY2Q_PACKEDBYTES];

    for (int t = 0; t < NTESTS; t++) {
        poly_uniform(&a, seed, nonce++);
        b = a;
        polyq_pack(rq, &a);
        polyq_unpack(&a, rq);
        for (int i = 0; i < N; ++i)
            ASSERT_EQ(a.coeffs[i], b.coeffs[i]);
        poly_uniform_eta(&a, seed, nonce++);
        b = a;
        polyeta_pack(reta, &a);
        polyeta_unpack(&a, reta);
        for (int i = 0; i < N; ++i)
            ASSERT_EQ(a.coeffs[i], b.coeffs[i]);
        for (int i = 0; i < N; i++) {
            a.coeffs[i] = rand();
        }
        poly_freeze2q(&a);
        b = a;
        poly2q_pack(r2q, &a);
        poly2q_unpack(&a, r2q);
        for (int i = 0; i < N; ++i)
            ASSERT_EQ(a.coeffs[i], b.coeffs[i]);
    }
}

void testPolyVecGen() {
    uint8_t seed[CRHBYTES];
    uint16_t nonce = 0;
    for (int t = 0; t < NTESTS; t++) {
        polyvecm b;
        randombytes(seed, sizeof(seed));
        polyvecmm_uniform_eta(&b, seed, nonce);
        nonce += M;
        for (int i = 0; i < L - 1; ++i) {
            for (int j = 0; j < N; ++j) {
                ASSERT_TRUE(-ETA <= b.vec[i].coeffs[j] &&
                            b.vec[i].coeffs[i] <= ETA);
            }
        }
    }
}

void testPolyMatGen() {
    uint8_t seed[SEEDBYTES];
    for (int t = 0; t < NTESTS; t++) {
        randombytes(seed, sizeof(seed));
        polyvecm a[K];
        polymat_expand(a, seed);
        if (PRINT_FLAG) {
            printf("A = ([");
            for (int j = 0; j < K; ++j) {
                for (int k = 0; k < M - 1; ++k) {
                    for (int l = 0; l < N; ++l) {
                        printf("%4d", a[j].vec[k].coeffs[l]);
                        if (l < N - 1)
                            printf(", ");
                        else if (k < M - 1)
                            printf("], [");
                        else if (j < K - 1)
                            printf("];\n     [");
                        else
                            printf("])\n");
                    }
                }
            }
        }
        for (int j = 0; j < K; ++j) {
            for (int k = 0; k < M - 1; ++k) {
                for (int l = 0; l < N; ++l) {
                    ASSERT_TRUE(a[j].vec[k].coeffs[l] >= 0 &&
                                a[j].vec[k].coeffs[l] < Q);
                }
            }
        }
    }
}

void testPolyMatInv() {
    uint8_t seed[CRHBYTES];
    uint16_t nonce = 0;
    randombytes(seed, sizeof(seed));
    for (int t = 0; t < NTESTS; t++) {
        polyveck a[K];
        polyveck b[K];
        polyveck c[K];
        polyveck d[K];
        do {
            for (int i = 0; i < K; ++i) {
                for (int j = 0; j < K; ++j) {
                    poly_uniform(&a[i].vec[j], seed, nonce++);
                }
                polyveck_ntt(&a[i]);
            }
            for (int i = 0; i < K; ++i)
                for (int j = 0; j < K; ++j)
                    c[i].vec[j] = a[i].vec[j];
        } while (!polymatkk_inverse(b, a));
        for (int i = 0; i < K; ++i) { // checks if the reduced matrix is Id
            for (int j = 0; j < K; ++j) {
                for (int k = 0; k < N; ++k) {
                    int val = i == j ? 1 : 0;
                    ASSERT_EQ((a[i].vec[j].coeffs[k] - val) % Q, 0);
                }
            }
        }
        poly tp;
        for (int i = 0; i < K; ++i) {
            for (int j = 0; j < K; ++j) {
                poly_pointwise_montgomery(&d[i].vec[j], &b[i].vec[0],
                                          &c[0].vec[j]);
                for (int k = 1; k < K; ++k) {
                    poly_pointwise_montgomery(&tp, &b[i].vec[k], &c[k].vec[j]);
                    poly_add(&d[i].vec[j], &d[i].vec[j], &tp);
                }
            }
        }
        for (int i = 0; i < K; ++i) {
            for (int j = 0; j < K; ++j) {
                poly_invntt_tomont(&d[i].vec[j]);
            }
        }
        for (int i = 0; i < K; ++i) { // checks if the multiplied mat is Id
            for (int j = 0; j < K; ++j) {
                for (int k = 0; k < N; ++k) {
                    int val = i == j && k == 0 ? 1 : 0;
                    ASSERT_EQ((d[i].vec[j].coeffs[k] - val) % Q, 0);
                }
            }
        }
    }
}

void testKeyPack() {
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t pk_copy[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk_copy[CRYPTO_SECRETKEYBYTES];
    uint8_t seed[SEEDBYTES];
    polyveck b;
    polyvecm A[K];
    polyvecm s;
    for (int t = 0; t < NTESTS; t++) {
        crypto_sign_keypair(pk, sk);
        memcpy(pk_copy, pk, CRYPTO_PUBLICKEYBYTES);
        memcpy(sk_copy, sk, CRYPTO_SECRETKEYBYTES);
        unpack_pk(seed, &b, pk_copy);
        pack_pk(pk_copy, seed, &b);
        ASSERT_TRUE(arrayEq(pk, pk_copy, CRYPTO_PUBLICKEYBYTES));
        unpack_sk(A, &s, sk_copy);
        pack_sk(sk_copy, pk_copy, &s);
        ASSERT_TRUE(arrayEq(sk, sk_copy, CRYPTO_SECRETKEYBYTES));
    }
}

void testKeyGen() { // tests if A * s' = 0 (mod q)
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    polyvecm A[K];
    polyvecm s;
    polyveck res;
    poly temp;
    for (int t = 0; t < NTESTS; t++) {
        crypto_sign_keypair(pk, sk);
        unpack_sk(A, &s, sk);
        for (int i = 0; i < N; ++i) {
            int val = i == 0 ? 1 : 0;
            ASSERT_EQ(val, s.vec[M - 1].coeffs[i]);
        }
        polymatkmm_naivemul(&res, A, &s);
        for (int i = 0; i < K; ++i) {
            poly_naivemul(&temp, &A[i].vec[M - 1], &s.vec[M - 1]);
            poly_add(&res.vec[i], &res.vec[i], &temp);
        }
        for (int i = 0; i < K; ++i) {
            for (int j = 0; j < N; ++j) {
                ASSERT_EQ((res.vec[i].coeffs[j]) % Q, 0); // qj
            }
        }
    }
}

void testSignPack() {
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sig[CRYPTO_BYTES];
    uint8_t sig_copy[CRYPTO_BYTES];
    uint8_t c[SEEDBYTES];
    uint8_t msg[SEEDBYTES];
    size_t siglen;
    polyvecm z;
    for (int t = 0; t < NTESTS; t++) {
        crypto_sign_keypair(pk, sk);
        randombytes(msg, SEEDBYTES);
        crypto_sign_signature(sig, &siglen, msg, SEEDBYTES, sk);
        memcpy(sig_copy, sig, CRYPTO_BYTES);
        unpack_sig(c, &z, sig_copy);
        pack_sig(sig_copy, c, &z);
        ASSERT_TRUE(arrayEq(sig, sig_copy, CRYPTO_BYTES));
    }
}

void testSignVerify() {
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sig[CRYPTO_BYTES];
    uint8_t msg[SEEDBYTES];
    size_t siglen;
    for (int t = 0; t < NTESTS; t++) {
        crypto_sign_keypair(pk, sk);
        randombytes(msg, SEEDBYTES);
        crypto_sign_signature(sig, &siglen, msg, SEEDBYTES, sk);
        ASSERT_EQ(crypto_sign_verify(sig, siglen, msg, SEEDBYTES, pk), 0);
    }
}

void testSignMismatch() {
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sig[CRYPTO_BYTES];
    uint8_t msg[SEEDBYTES];
    size_t siglen;
    for (int t = 0; t < NTESTS; t++) {
        crypto_sign_keypair(pk, sk);
        randombytes(msg, SEEDBYTES);
        crypto_sign_signature(sig, &siglen, msg, SEEDBYTES, sk);
        mutateSignature(sig);
        ASSERT_EQ(crypto_sign_verify(sig, siglen, msg, SEEDBYTES, pk), -1);
    }
}

TEST(POLY, GEN) { testPolyGen(); }

TEST(POLY, NTT) { testPolyNTT(); }

TEST(POLY, MULT) { testPolyMult(); }

TEST(POLY, INV) { testPolyInv(); }

TEST(POLY, CRT) { testPolyCRT(); }

TEST(POLY, PACK) { testPolyPack(); }

TEST(POLYVEC, GEN) { testPolyVecGen(); }

TEST(POLYMAT, GEN) { testPolyMatGen(); }

TEST(POLYMAT, INV) { testPolyMatInv(); }

TEST(KEY, PACK) { testKeyPack(); }

TEST(KEY, GEN) { testKeyGen(); }

TEST(SIGN, PACK) { testSignPack(); }

TEST(SIGN, VERIFY) { testSignVerify(); }

TEST(SIGN, MISMATCH) { testSignMismatch(); }