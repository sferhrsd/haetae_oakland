#include "../include/fft.h"
#include <stdio.h>
#include <stdint.h>

int main(void)
{
  int i;
  complex_fp32_16 input[FFT_N] = {0};
  input[128].real = 1<<16;
  input[384].real = 1<<16;

  fft(input);

  double p16 = 1UL<<16;
  for (i = 0; i < FFT_N; i++)
  {
    printf("%f, %f\n", input[i].real / p16, input[i].imag / p16);
  }

  return 0;
}
