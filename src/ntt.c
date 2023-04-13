#include "ntt.h"
#include "params.h"
#include "reduce.h"
#include <stdint.h>

// TODO: different zetas / f according to Q

static const int32_t zetas[N] = {
     0,  26964, -16505,  22229,  30746,  20243,  19064, -31218, 
  9395, -30985,  22859,  -8851,  32144,  13744,  21408,  17599, 
-16039, -22946,   6241, -19553,  10681,  22935,  22431, -29104, 
 28147, -27527, -29133, -20035,  20143, -11361,  30820,  25252, 
-22562,  -6789, -10049,   9383,  16304, -12296,  16446,  18239, 
 -1296, -19725, -32076,  11782, -17941,  29643,  -8577,   7893, 
-21464, -19646, -15130,  -2391,  30608, -23970, -16608,  19616, 
 -7941,  26533, -19129,  27690,   7597, -11459,  10615,  -9430, 
 11591,   7814,  12697,  32114,  -3761,  -9604,  19813,  20353, 
 17456, -16267, -19555,    598, -29942,   4538,    835,  15546, 
  3970, -27685,   1488,   8311, -12442,  31352, -17631,   1806, 
 -5342,   9790,  29068,  16507, -29051,  22131,   6759,  15510, 
-14941,  28710,   1160, -31327,  24985,  11261, -10623, -27727, 
 21502,  18731, -16186,  -4127, -18832,  12050, -14501,   7929, 
 29563, -31064,   5913,   5322, -16405,   2844,  29439,   5876, 
 -9522, -18586,  -9874,  23844,  30362, -21442,   9560,  17671, 
-27989,   3350,    787, -13857,   1657, -21224,  -7374,  -9190, 
  2464,  25555,  -3529, -28772,  16588, -15739,  23475,  13666, 
  5764,  30980,  13633,  -7401, -30317,  28847,   7682, -11808, 
 -8796,  14864, -24162, -19194,    689,  -1311, -31332, -16319, 
  1025,  10971, -23016,  -2648, -21900, -12543, -25921,  28254, 
 28521, -16160,  12380, -12882, -30332, -16630,  23439,   7742, 
 17182,  17494,   5920,  13642,   7382, -18166,  21422, -30274, 
-28190,  13283, -20316,  -9939,  10672,  21454,   6080, -17374, 
-29735, -25912, -10170,   3808,  10639, -26985, -10865,  25636, 
 17261, -26851,  -8253,  -3304,  18282,  -2202, -31368, -22243, 
 13882,  12069, -11242,  -7729, -10226,   1761, -27298,  -4800, 
-17737, -22805,  -3528,     65,  10770,   8908, -23751,  26934, 
 21921, -27010, -21944,   8889,  -1035,  23224,  -9488,  -5823, 
  -994, -20206,   7655, -16251, -22820, -27740,  15822,  23078, 
 13803,  -8099,   2931,   9217, -21126, -14203,  25492, -12831, 
  7947,  17463, -12979,  29003,  31612,  26554,   8241, -20175}; // q = 64513

/*************************************************
 * Name:        ntt
 *
 * Description: Forward NTT, in-place. No modular reduction is performed after
 *              additions or subtractions. Output vector is in bitreversed
 *order.
 *
 * Arguments:   - uint32_t p[N]: input/output coefficient array
 **************************************************/
void ntt(int32_t a[N]) {
    unsigned int len, start, j, k;
    int32_t zeta, t;

    k = 0;
    for (len = 128; len > 0; len >>= 1) {
        for (start = 0; start < N; start = j + len) {
            zeta = zetas[++k];
            for (j = start; j < start + len; ++j) {
                t = montgomery_reduce((int64_t)zeta * a[j + len]);
                a[j + len] = a[j] - t;
                a[j] = a[j] + t;
            }
        }
    }
}

/*************************************************
 * Name:        invntt_tomont
 *
 * Description: Inverse NTT and multiplication by Montgomery factor 2^32.
 *              In-place. No modular reductions after additions or
 *              subtractions; input coefficients need to be smaller than
 *              Q in absolute value. Output coefficient are smaller than Q in
 *              absolute value.
 *
 * Arguments:   - uint32_t p[N]: input/output coefficient array
 **************************************************/
void invntt_tomont(int32_t a[N]) {
    unsigned int start, len, j, k;
    int32_t t, zeta;
    const int32_t f = -29720; // mont^2/256

    k = 256;
    for (len = 1; len < N; len <<= 1) {
        for (start = 0; start < N; start = j + len) {
            zeta = -zetas[--k];
            for (j = start; j < start + len; ++j) {
                t = a[j];
                a[j] = t + a[j + len];
                a[j + len] = t - a[j + len];
                a[j + len] = montgomery_reduce((int64_t)zeta * a[j + len]);
            }
        }
    }

    for (j = 0; j < N; ++j) {
        a[j] = montgomery_reduce((int64_t)f * a[j]);
    }
}
