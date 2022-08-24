#include <cassert>
#include <cmath>
#include <cstdlib>
#include "tfhe_core.h"
#include "polynomials_arithmetic.h"
#include "lagrangehalfc_arithmetic.h"

#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#include <iostream>
#include "2.h"

//2.h = 1024 coef  3 =2048 4=4096



using namespace std;

#define mod_add(a, b, modulus) ((modulus - a) > b) ? (a + b) : (a + b - modulus)
#define mod_sub(a, b, modulus) (a >= b) ? (a - b) : (modulus - b + a)





uint64_t mod_mul_p(uint64_t x, uint64_t y, uint64_t p) {




    uint32_t x0 = (uint32_t)x;
    uint32_t x1 = (uint32_t)(x >> 32);
    uint32_t y0 = (uint32_t)y;
    uint32_t y1 = (uint32_t)(y >> 32);




    uint64_t x0y0 = (uint64_t)x0 * (uint64_t)y0;
    uint64_t x0y1 = (uint64_t)x0 * (uint64_t)y1;
    uint64_t x1y0 = (uint64_t)x1 * (uint64_t)y0;
    uint64_t x1y1 = (uint64_t)x1 * (uint64_t)y1;



    uint32_t d = (uint32_t)x0y0;
    uint64_t pp1 = (x0y0 >> 32) + (uint32_t)(x1y0)+(uint32_t)(x0y1);
    uint32_t c = (uint32_t)pp1;
    uint64_t pp2 = (x1y0 >> 32) + (x0y1 >> 32) + (uint32_t)(x1y1);
    uint64_t pp3 = (pp1 >> 32) + (uint32_t)(pp2);

    uint32_t a = (pp2 >> 32) + (x1y1 >> 32);

    uint64_t bpc = (uint32_t)pp3 + (uint64_t)c; //1fffffffe max




    bpc = ((bpc + (bpc >> 32)) << 32) - (bpc >> 32); //1fffffffe +

    uint64_t mines = ((uint64_t)a + ((uint64_t)(uint32_t)pp3));
    uint64_t plus = bpc + (uint64_t)d;

    if (plus >= mines)
        return (plus - mines);
    else
        return p - mines + plus;

    /* result = result +p;
     if(result>p)
         result -=p;
    */






}







//////////////////////////////////////////////////////////////////////////






void bitrev_shuffle_p(uint64_t x[], int N) {

    int j = 0;
    for (int i = 1; i < N; i++) {
        int b = N >> 1;
        while (j >= b) {
            j -= b;
            b >>= 1;
        }
        j += b;
        if (j > i) {
            uint64_t temp = x[j];
            x[j] = x[i];
            x[i] = temp;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

void ntt_CT_p(uint64_t x[], uint64_t modulus, int N) {
    uint64_t b = 0;
    //uint64_t c=0;

    bitrev_shuffle_p(x, N);

    int wbarr = 0;
    for (int trans_size = 2; trans_size <= N; trans_size = trans_size * 2) {
        uint64_t wb = 1;
        for (int t = 0; t < (trans_size >> 1); t++) {
            for (int trans = 0; trans < (N / trans_size); trans++) {
                int i = trans * trans_size + t;

                int j = i + (trans_size >> 1);
                uint64_t a = x[i];
                if (wb == 1) {
                    b = x[j];
                }
                else {

                    //c=mod_mul_p2(x[j], wb, modulus);

                    b = mod_mul_p(x[j], wb, modulus);




                }
                x[i] = mod_add(a, b, modulus);
                x[j] = mod_sub(a, b, modulus);


            }
            wb = wb_normal_2[wbarr];
            wbarr++;

        }
    }

}




void intt_CT_p(uint64_t x[], uint64_t modulus, uint64_t y[], int N) {
    uint64_t b = 0;
    int inv = 0;

    bitrev_shuffle_p(x, N);


    for (int trans_size = 2; trans_size <= N; trans_size = trans_size * 2) {
        uint64_t wb = 1;


        for (int t = 0; t < (trans_size >> 1); t++) {
            for (int trans = 0; trans < (N / trans_size); trans++) {
                int i = trans * trans_size + t;
                int j = i + (trans_size >> 1);
                uint64_t a = x[i];


                if (wb == 1)
                    b = x[j];
                else
                    b = mod_mul_p(x[j], wb, modulus);

                x[i] = mod_add(a, b, modulus);
                x[j] = mod_sub(a, b, modulus);

            }
            wb = wb_inverse_2[inv];
            inv++;
        }
    }



    int i = 0;

    for (i = 0; i < N; i++) {
        y[i] = mod_mul_p(x[i], scale_2, modulus);    //scale_2*phi inversev, modulus  (phi inverse sclaed)
        uint64_t phi = phi_inverse_2[i];

        y[i] = mod_mul_p(y[i], phi, modulus);
    }

}





//allocate memory space for a LagrangeHalfCPolynomial
EXPORT LagrangeHalfCPolynomial* alloc_LagrangeHalfCPolynomial() {
    return (LagrangeHalfCPolynomial*)malloc(sizeof(LagrangeHalfCPolynomial));
}
EXPORT LagrangeHalfCPolynomial* alloc_LagrangeHalfCPolynomial_array(int32_t nbelts) {
    return (LagrangeHalfCPolynomial*)malloc(nbelts * sizeof(LagrangeHalfCPolynomial));
}

//free memory space for a LweKey
EXPORT void free_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* ptr) {
    free(ptr);
}
EXPORT void free_LagrangeHalfCPolynomial_array(int32_t nbelts, LagrangeHalfCPolynomial* ptr) {
    free(ptr);
}

//allocates and initialize the LagrangeHalfCPolynomial structure
//(equivalent of the C++ new)
EXPORT LagrangeHalfCPolynomial* new_LagrangeHalfCPolynomial(const int32_t N) {
    LagrangeHalfCPolynomial* obj = alloc_LagrangeHalfCPolynomial();
    init_LagrangeHalfCPolynomial(obj, N);
    return obj;
}
EXPORT LagrangeHalfCPolynomial* new_LagrangeHalfCPolynomial_array(int32_t nbelts, const int32_t N) {
    LagrangeHalfCPolynomial* obj = alloc_LagrangeHalfCPolynomial_array(nbelts);
    init_LagrangeHalfCPolynomial_array(nbelts, obj, N);
    return obj;
}

//destroys and frees the LagrangeHalfCPolynomial structure
//(equivalent of the C++ delete)
EXPORT void delete_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* obj) {
    destroy_LagrangeHalfCPolynomial(obj);
    free_LagrangeHalfCPolynomial(obj);
}
EXPORT void delete_LagrangeHalfCPolynomial_array(int32_t nbelts, LagrangeHalfCPolynomial* obj) {
    destroy_LagrangeHalfCPolynomial_array(nbelts, obj);
    free_LagrangeHalfCPolynomial_array(nbelts, obj);
}

/** multiplication via direct FFT (it must know the implem of LagrangeHalfCPolynomial because of the tmp+1 notation */
EXPORT void torusPolynomialMultFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2) {
    const int32_t N =1024; // N=1024

      /*static int z = 0;
      printf("%d invocation\n",z++); */

    uint64_t  modulus = 0xffffffff00000001UL; //


    //uint64_t w = 4440654710286119610; // N=1024
    //printf("w      : %" PRIu64 "\n", w);


   // uint64_t inv_w = 12337821426711963180U; // N=1024

  //printf("%d invocation torusPolynomialMultFFT\n",zz++);

  //LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3,N);
  //IntPolynomial_ifft(tmp+0,poly1);
  //TorusPolynomial_ifft(tmp+1,poly2);
  //LagrangeHalfCPolynomialMul(tmp+2,tmp+0,tmp+1);
  //TorusPolynomial_fft(result, tmp+2); 
  //delete_LagrangeHalfCPolynomial_array(3,tmp);

    uint64_t a1[N], a2[N], intt_mul_res[N], ntt_res[N]/*, intt_a1_res[N], intt_a2_res[N]*/;

    for (int i = 0; i < N; i++) {
        // a1[i] = (uint64_t)poly1->coefs[i] ;
        if (poly1->coefs[i] < 0) {
            //uint64_t phi = 9630642590298920985U; // 8816101479115663336U;
            a1[i] = poly1->coefs[i] + modulus;
            uint64_t phi = phi_normal_2[i];
            a1[i] = mod_mul_p(a1[i], phi, modulus);
        }
        else {
            //uint64_t phi = 9630642590298920985U; // 8816101479115663336U;
            a1[i] = (uint64_t)poly1->coefs[i];
            uint64_t phi = phi_normal_2[i];
            a1[i] = mod_mul_p(a1[i], phi, modulus);
        }

        //printf(" ntt poly1->coefs[i] %d   : %d\n", i,poly1->coefs[i]);
        //printf(" ntt a1 %d   : %" PRIu64 "\n", i,a1[i]);
        //a2[i] = (uint64_t)poly2->coefsT[i] ;
        if (poly2->coefsT[i] < 0) {
            //uint64_t phi = 9630642590298920985U; // 8816101479115663336U;
            a2[i] = poly2->coefsT[i] + modulus;
            uint64_t phi = phi_normal_2[i];//mod_exp(phi, i, modulus);
            a2[i] = mod_mul_p(a2[i], phi, modulus);
        }
        else {
            //uint64_t phi = 9630642590298920985U; // 8816101479115663336U;
            a2[i] = (uint64_t)poly2->coefsT[i];
            uint64_t phi = phi_normal_2[i];
            a2[i] = mod_mul_p(a2[i], phi, modulus);
        }

        //printf(" ntt a2 %d   : %" PRIu64 "\n", i,a2[i]);
        /*intt_a1_res[i] = 0;
        intt_a2_res[i] = 0;
        ntt_res[i] = 0;*/
    }



    ntt_CT_p(a1, modulus, N);


    ntt_CT_p(a2, modulus, N);


    for (int32_t i = 0; i < N; i++)
        intt_mul_res[i] = mod_mul_p(a1[i], a2[i], modulus);



    intt_CT_p(intt_mul_res, modulus, ntt_res, N);


    uint64_t med = modulus / 2;
    for (int i = 0; i < N; i++) {
        result->coefsT[i] = (Torus32)((ntt_res[i] & 0xffffffff) - (ntt_res[i] > med)); // worked!


    }


}
/** multiplication via direct FFT (it must know the implem of LagrangeHalfCPolynomial because of the tmp+1 notation */
//EXPORT void torusPolynomialMultFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2) {
EXPORT void torusPolynomialAddMulRFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2) {
    const int32_t N =1024; // N=1024

      /*static int z = 0;
      printf("%d invocation\n",z++); */

    uint64_t  modulus = 0xffffffff00000001UL; //


    //uint64_t w = 4440654710286119610; // N=1024
    //printf("w      : %" PRIu64 "\n", w);


   // uint64_t inv_w = 12337821426711963180U; // N=1024

  //printf("%d invocation torusPolynomialMultFFT\n",zz++);

  //LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3,N);
  //IntPolynomial_ifft(tmp+0,poly1);
  //TorusPolynomial_ifft(tmp+1,poly2);
  //LagrangeHalfCPolynomialMul(tmp+2,tmp+0,tmp+1);
  //TorusPolynomial_fft(result, tmp+2);
  //delete_LagrangeHalfCPolynomial_array(3,tmp);

    uint64_t a1[N], a2[N], intt_mul_res[N], ntt_res[N]/*, intt_a1_res[N], intt_a2_res[N]*/;

    for (int i = 0; i < N; i++) {
        // a1[i] = (uint64_t)poly1->coefs[i] ;
        if (poly1->coefs[i] < 0) {
            //uint64_t phi = 9630642590298920985U; // 8816101479115663336U;
            a1[i] = poly1->coefs[i] + modulus;
            uint64_t phi = phi_normal_2[i];
            a1[i] = mod_mul_p(a1[i], phi, modulus);
        }
        else {
            //uint64_t phi = 9630642590298920985U; // 8816101479115663336U;
            a1[i] = (uint64_t)poly1->coefs[i];
            uint64_t phi = phi_normal_2[i];
            a1[i] = mod_mul_p(a1[i], phi, modulus);
        }

        //printf(" ntt poly1->coefs[i] %d   : %d\n", i,poly1->coefs[i]);
        //printf(" ntt a1 %d   : %" PRIu64 "\n", i,a1[i]);
        //a2[i] = (uint64_t)poly2->coefsT[i] ;
        if (poly2->coefsT[i] < 0) {
            //uint64_t phi = 9630642590298920985U; // 8816101479115663336U;
            a2[i] = poly2->coefsT[i] + modulus;
            uint64_t phi = phi_normal_2[i];//mod_exp(phi, i, modulus);
            a2[i] = mod_mul_p(a2[i], phi, modulus);
        }
        else {
            //uint64_t phi = 9630642590298920985U; // 8816101479115663336U;
            a2[i] = (uint64_t)poly2->coefsT[i];
            uint64_t phi = phi_normal_2[i];
            a2[i] = mod_mul_p(a2[i], phi, modulus);
        }

        //printf(" ntt a2 %d   : %" PRIu64 "\n", i,a2[i]);
        /*intt_a1_res[i] = 0;
        intt_a2_res[i] = 0;
        ntt_res[i] = 0;*/
    }



    ntt_CT_p(a1, modulus, N);


    ntt_CT_p(a2, modulus, N);


    for (int32_t i = 0; i < N; i++)
        intt_mul_res[i] = mod_mul_p(a1[i], a2[i], modulus);



    intt_CT_p(intt_mul_res, modulus, ntt_res, N);

    TorusPolynomial* tmpr = new_TorusPolynomial(N);
    uint64_t med = modulus / 2;
    for (int i = 0; i < N; i++) {
        tmpr->coefsT[i] = (Torus32)((ntt_res[i] & 0xffffffff) - (ntt_res[i] > med)); // worked!
    }

    torusPolynomialAddTo(result, tmpr);

}

EXPORT void torusPolynomialSubMulRFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2) {

    const int32_t N = 1024; // N=1024

      /*static int z = 0;
      printf("%d invocation\n",z++); */

    uint64_t  modulus = 0xffffffff00000001UL; //


    //uint64_t w = 4440654710286119610; // N=1024
    //printf("w      : %" PRIu64 "\n", w);


   // uint64_t inv_w = 12337821426711963180U; // N=1024

  //printf("%d invocation torusPolynomialMultFFT\n",zz++);

  //LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3,N);
  //IntPolynomial_ifft(tmp+0,poly1);
  //TorusPolynomial_ifft(tmp+1,poly2);
  //LagrangeHalfCPolynomialMul(tmp+2,tmp+0,tmp+1);
  //TorusPolynomial_fft(result, tmp+2);
  //delete_LagrangeHalfCPolynomial_array(3,tmp);

    uint64_t a1[N], a2[N], intt_mul_res[N], ntt_res[N]/*, intt_a1_res[N], intt_a2_res[N]*/;

    for (int i = 0; i < N; i++) {
        // a1[i] = (uint64_t)poly1->coefs[i] ;
        if (poly1->coefs[i] < 0) {
            //uint64_t phi = 9630642590298920985U; // 8816101479115663336U;
            a1[i] = poly1->coefs[i] + modulus;
            uint64_t phi = phi_normal_2[i];
            a1[i] = mod_mul_p(a1[i], phi, modulus);
        }
        else {
            //uint64_t phi = 9630642590298920985U; // 8816101479115663336U;
            a1[i] = (uint64_t)poly1->coefs[i];
            uint64_t phi = phi_normal_2[i];
            a1[i] = mod_mul_p(a1[i], phi, modulus);
        }

        //printf(" ntt poly1->coefs[i] %d   : %d\n", i,poly1->coefs[i]);
        //printf(" ntt a1 %d   : %" PRIu64 "\n", i,a1[i]);
        //a2[i] = (uint64_t)poly2->coefsT[i] ;
        if (poly2->coefsT[i] < 0) {
            //uint64_t phi = 9630642590298920985U; // 8816101479115663336U;
            a2[i] = poly2->coefsT[i] + modulus;
            uint64_t phi = phi_normal_2[i];//mod_exp(phi, i, modulus);
            a2[i] = mod_mul_p(a2[i], phi, modulus);
        }
        else {
            //uint64_t phi = 9630642590298920985U; // 8816101479115663336U;
            a2[i] = (uint64_t)poly2->coefsT[i];
            uint64_t phi = phi_normal_2[i];
            a2[i] = mod_mul_p(a2[i], phi, modulus);
        }

        //printf(" ntt a2 %d   : %" PRIu64 "\n", i,a2[i]);
        /*intt_a1_res[i] = 0;
        intt_a2_res[i] = 0;
        ntt_res[i] = 0;*/
    }



    ntt_CT_p(a1, modulus, N);


    ntt_CT_p(a2, modulus, N);


    for (int32_t i = 0; i < N; i++)
        intt_mul_res[i] = mod_mul_p(a1[i], a2[i], modulus);



    intt_CT_p(intt_mul_res, modulus, ntt_res, N);

    TorusPolynomial* tmpr = new_TorusPolynomial(N);
    uint64_t med = modulus / 2;
    for (int i = 0; i < N; i++) {
        tmpr->coefsT[i] = (Torus32)((ntt_res[i] & 0xffffffff) - (ntt_res[i] > med)); // worked!
    }


    torusPolynomialSubTo(result, tmpr);

}
