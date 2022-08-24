#include<stdint.h>
#include<stdio.h>
#include "const_1024.h"



 uint64_t modular_mult_c(uint64_t x, uint64_t y, uint64_t p);
 void ntt_c(uint64_t x[], uint64_t modulus, int N) ;
 void intt_c(uint64_t x[], uint64_t modulus, uint64_t y[], int N);
 int polymuladd(uint64_t Instream[8192],int32_t outStream[1024], int batu);
