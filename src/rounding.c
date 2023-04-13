#include "rounding.h"
#include "params.h"
#include <stdint.h>

/*************************************************
 * Name:        decompose
 *
 * Description: For finite field element r, compute high and lowbits 
 *              hb, lb such that r = hb * b + lb with -b/4 < lb <= b/4.
 *
 * Arguments:   - int32_t r: input element
 *              - int32_t *lowbits: pointer to output element lb
 *              - int32_t *highbits: pointer to output element hb
 *
 * Returns a1.
 **************************************************/
void decompose(int32_t *highbits, int32_t *lowbits, const int32_t r) {
    int32_t lb, center;
    uint32_t alpha_mask = BETA-1;

    lb = r & alpha_mask;
    center = (HALF_BETA - (lb+1)) >>31; // if lb >= HALF_ALPHA
    lb -= BETA & center;
    *lowbits = lb;

    *highbits = (r + HALF_BETA) >> LOG_BETA;
}

void decompose_hint(int32_t *highbits, int32_t *lowbits, const int32_t r) {
    int32_t lb, hb, center, edgecase;
    uint32_t alpha_mask = ALPHA-1;

    //TODO: remove lowbits?
    lb = r & alpha_mask;
    center = (HALF_ALPHA - (lb+1)) >>31; // if lb >= HALF_ALPHA
    lb -= ALPHA & center;

    hb = (r + HALF_ALPHA) >> LOG_ALPHA;

    edgecase = ((DQ-2)/ALPHA - (hb+1)) >>31; // if hb == (DQ-2)/ALPHA
    lb -= 2 & edgecase;    // lb = lb - 2
    hb -= (DQ-2)/ALPHA & edgecase; // hb = 0

    *lowbits = lb;
    *highbits = hb;
}

/*************************************************
* Name:        power2round
*
* Description: For finite field element a, compute a0, a1 such that
*              a mod^+ Q = a1*2^D + a0 with -2^{D-1} < a0 <= 2^{D-1}.
*              Assumes a to be standard representative.
*
* Arguments:   - int32_t a: input element
*              - int32_t *a0: pointer to output element a0
*
* Returns a1.
**************************************************/
int32_t power2round(int32_t *a0, const int32_t a)  {
#if D > 1
#error "Only implemented for D = 1"
#endif
  *a0 = a&1;
  *a0 -= ((a>>1)&*a0) << 1;
  return (a-*a0) >> 1;
}
