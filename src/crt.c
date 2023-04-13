#include "crt.h"
#include "params.h"
#include <stdint.h>

int32_t fromcrt(int32_t xq, int32_t x2)
{
  return xq + (Q & -((xq^x2)&1));
}

int32_t fromcrt0(int32_t xq)
{
  return xq + (Q & -(xq&1));
}

