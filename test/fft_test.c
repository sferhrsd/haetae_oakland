#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <inttypes.h>

#include "fftw3.h"
#include "fips202.h"
#include "params.h"
#include "polyvec.h"
#include "randombytes.h"

#define N0 (1 << 5)
#define num_test (10000)

int fft_basic_test(void);
int fft_sum(void);
int fft_execute_test(void);
int singular_value_test(void);

void print_complex(fftw_complex *input, size_t input_len) {
    for (size_t i = 0; i < input_len; ++i) {
        printf("(%lf, %lf)\n", creal(input[i]), cimag(input[i]));
    }
}

int main(void) {
    printf("HAETAE mode = %d with gamma = %d\n", HAETAE_MODE, GAMMA);

    size_t print_count = 1;
    uint8_t entropy[48] = {0};
    randombytes_init(entropy, NULL, 256);

    int fail_count = 0;
    for (int i = 0; i < num_test; ++i) {
        // if (fft_basic_test()) {
        //     printf("fft test failure\n");
        // }

        // if (fft_sum()) {
        //     printf("fft sum failure\n");
        // }

        // if (fft_execute_test()) {
        //     printf("fft execute failure\n");
        // }

        if (singular_value_test()) {
            printf("At %d, singular value test failure\n", i);
            ++fail_count;
        }

        if (!(i % (num_test / 10))) {
            printf("...%lu%%", print_count * 10);
            fflush(stdout);
            ++print_count;
        }
    }

    if (fail_count) {
        printf("Fail %d cases of %d total\n", fail_count, num_test);
    }

    return 0;
}

int fft_basic_test(void) {
    int i = 0;
    int input[N0] = {0};
    int64_t inint[N0], resint[N0];
    double *in, *res;
    fftw_complex *out;
    fftw_plan p;

    in = (double *)fftw_malloc(sizeof(double) * N0);
    res = (double *)fftw_malloc(sizeof(double) * N0);
    out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N0);
    out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N0);

    randombytes((uint8_t *)input, sizeof(int) * N0);

    for (i = 0; i < N0; ++i) {
        in[i] = input[i] % 100;
        inint[i] = in[i];
    }

    // foward fft
    p = fftw_plan_dft_r2c_1d(N0, in, out, FFTW_ESTIMATE);
    fftw_execute(p); /* repeat as needed */
    // print_complex(out, N0);

    // backward fft
    p = fftw_plan_dft_c2r_1d(N0, out, res, FFTW_ESTIMATE);
    fftw_execute(p); /* repeat as needed */

    // compare
    for (i = 0; i < N0; ++i) {
        res[i] /= N0;
        resint[i] = lround(res[i]);
        if (inint[i] != resint[i]) {
            printf("Diff at %d with ", i);
            printf("inint[%d]=%ld and resint[%d]=%ld and ", i, inint[i], i,
                   resint[i]);
            printf("res[%d]=%lf\n", i, res[i]);
            break;
        }
    }

    // free
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(res);
    fftw_free(out);

    if (i == N0) {
        printf("FFT basic test good\n");
        return 0;
    } else {
        return 1;
    }
}

int fft_sum(void) {
    int i = 0;
    int32_t answer[N0], resint[N0];
    double *in1, *in2, *res;
    fftw_complex *out1, *out2;
    fftw_plan p;

    in1 = (double *)malloc(sizeof(double) * N0);
    in2 = (double *)malloc(sizeof(double) * N0);
    res = (double *)malloc(sizeof(double) * N0);
    out1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N0);
    out2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N0);

    for (int i = 0; i < N0; ++i) {
        in1[i] = i;
        in2[i] = i;
        answer[i] = 2 * i;
    }

    // forward fft
    p = fftw_plan_dft_r2c_1d(N0, in1, out1, FFTW_ESTIMATE);
    fftw_execute(p); /* repeat as needed */
    p = fftw_plan_dft_r2c_1d(N0, in2, out2, FFTW_ESTIMATE);
    fftw_execute(p); /* repeat as needed */

    // fft sum
    for (i = 0; i < N0; ++i) {
        out2[i] += out1[i];
    }

    // backward fft
    p = fftw_plan_dft_c2r_1d(N0, out2, res, FFTW_ESTIMATE);
    fftw_execute(p); /* repeat as needed */

    for (i = 0; i < N0; ++i) {
        res[i] /= N0;
        resint[i] = lround(res[i]);
    }

    for (i = 0; i < N0; ++i) {
        if (answer[i] != resint[i])
            printf("Diff at %d with answer[%d]=%d and resint[%d]=%d\n", i, i,
                   answer[i], i, resint[i]);
    }

    fftw_destroy_plan(p);
    free(in1);
    free(in2);
    free(res);
    fftw_free(out1);
    fftw_free(out2);

    if (i == N0) {
        printf("FFT sum good\n");
        return 0;
    } else
        return 1;
}

int fft_execute_test(void) {
    int i = 0;
    int64_t inint[N0], resint[N0];
    double *in, *res;
    fftw_complex *out;
    fftw_plan p;

    in = (double *)fftw_malloc(sizeof(double) * N0);
    res = (double *)fftw_malloc(sizeof(double) * N0);
    out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N0);

    for (i = 0; i < N0; ++i) {
        in[i] = i;
        inint[i] = in[i];
    }

    // foward fft
    p = fftw_plan_dft_r2c_1d(N0, in, out, FFTW_ESTIMATE);
    fftw_execute(p); /* repeat as needed */
    // print_complex(out, N0);

    int ran[N0] = {0};
    randombytes((uint8_t *)ran, N0 * sizeof(int));

    for (i = 0; i < N0; ++i) {
        in[i] = i + ran[i] % 10;
        inint[i] = in[i];
    }
    fftw_execute(p); /* repeat as needed */
    // print_complex(out, N0);

    // backward fft
    p = fftw_plan_dft_c2r_1d(N0, out, res, FFTW_ESTIMATE);
    fftw_execute(p); /* repeat as needed */

    // compare
    for (i = 0; i < N0; ++i) {
        res[i] /= N0;
        resint[i] = lround(res[i]);
        if (inint[i] != resint[i]) {
            printf("Diff at %d with in[%d]=%ld and res[%d]=%ld\n", i, i,
                   inint[i], i, resint[i]);
            break;
        }
    }

    // free
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(res);
    fftw_free(out);

    if (i == N0) {
        printf("FFT basic test good\n");
        return 0;
    } else {
        return 1;
    }
}

int singular_value_test(void) {
    uint8_t seedbuf[CRHBYTES] = {0};
    const uint8_t *sigma;
    polyvecm s1;
    polyveck s2;
    int count = 0;

reject:
    if (count == N0)
        return 1;

    // Get entropy \rho
    randombytes(seedbuf, CRHBYTES);

    // Sample seeds with entropy \rho
    shake256(seedbuf, CRHBYTES, seedbuf, CRHBYTES);
    sigma = seedbuf;

    // Sample secret vectors s1 and s2
    polyvecmk_uniform_eta(&s1, &s2, sigma);
    int64_t singular_value = polyvecmk_sqsing_value(&s1, &s2);
    printf("%" PRId64 "\n", singular_value);
    if ((singular_value > GAMMA * GAMMA)) {
        // printf("singular value = %lu\n", singular_value);
        ++count;
        goto reject;
    }

    return 0;
}
