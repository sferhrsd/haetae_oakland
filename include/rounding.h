#ifndef HAETAE_ROUNDING_H
#define HAETAE_ROUNDING_H

#include "params.h"
#include <stdint.h>

#define decompose HAETAE_NAMESPACE(decompose)
void decompose(int32_t *r2, int32_t *r1, const int32_t r);

#define power2round HAETAE_NAMESPACE(power2round)
int32_t power2round(int32_t *a0, const int32_t a);

#define decompose_hint HAETAE_NAMESPACE(decompose_hint)
void decompose_hint(int32_t *r2, int32_t *r1, const int32_t r);

#endif
