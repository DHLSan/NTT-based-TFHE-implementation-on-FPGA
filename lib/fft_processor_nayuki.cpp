#include <complex>
#include <polynomials.h>
#include "lagrangehalfc_impl.h"
#include "fft.h"
#include <cassert>
#include <cmath>
#include<iostream>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "const_1024.h"

#define mod_add(a, b, modulus) ((modulus - a) > b) ? (a + b) : (a + b - modulus)
#define mod_sub(a, b, modulus) (a >= b) ? (a - b) : (modulus - b + a)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
FFT_Processor_nayuki::FFT_Processor_nayuki(const int32_t N) : _2N(2 * N), N(N), Ns2(N / 2) {

}


void FFT_Processor_nayuki::check_alternate_real() {

}

void FFT_Processor_nayuki::execute_reverse_int(uint64_t* res, const IntPolynomial* a) {

    uint64_t modulus = 0xffffffff00000001UL;
    int N = a->N;

    for (int i = 0; i < N; i++) {
        if (a->coefs[i] < 0) {
            res[i] = uint64_t(a->coefs[i]) + modulus;
            uint64_t phi = phi_normal[i];
            res[i] = mod_mul(res[i], phi, modulus);
        }
        else {

            res[i] = (uint64_t)a->coefs[i];
            uint64_t phi = phi_normal[i];
            res[i] = mod_mul(res[i], phi, modulus);
        }
    }
    ntt_CT(res, modulus, N);




}

void FFT_Processor_nayuki::execute_reverse_torus32(uint64_t* res, const TorusPolynomial* a) {
    uint64_t modulus = 0xffffffff00000001UL;
    int N = a->N;

    for (int i = 0; i < N; i++) {
        if (a->coefsT[i] < 0) {
            res[i] = uint64_t(a->coefsT[i]) + modulus;
            uint64_t phi = phi_normal[i];
            res[i] = mod_mul(res[i], phi, modulus);
        }
        else {

            res[i] = (uint64_t)a->coefsT[i];
            uint64_t phi = phi_normal[i];
            res[i] = mod_mul(res[i], phi, modulus);
        }
    }

    ntt_CT(res, modulus, N);
}

void FFT_Processor_nayuki::execute_direct_torus32(TorusPolynomial* res, uint64_t* a) {
    uint64_t modulus = 0xffffffff00000001UL;
    int N = res->N;
    uint64_t intt_mul_res[1024];
    intt_CT(a, modulus, intt_mul_res, N);

    TorusPolynomial* tmpr = new_TorusPolynomial(N);

    uint64_t med = modulus / 2;
    for (int i = 0; i < N; i++) {
        tmpr->coefsT[i] =
            (Torus32)((intt_mul_res[i] & 0xffffffff) - (intt_mul_res[i] > med)); // worked!
    }

    for (int i = 0; i < 1024; i++) {
        res->coefsT[i] = tmpr->coefsT[i];
    }

}

FFT_Processor_nayuki::~FFT_Processor_nayuki() {

}

thread_local FFT_Processor_nayuki fp1024_nayuki(1024);

/**
 * FFT functions
 */
EXPORT void IntPolynomial_ifft(LagrangeHalfCPolynomial* result, const IntPolynomial* p) {
    LagrangeHalfCPolynomial_IMPL* r = (LagrangeHalfCPolynomial_IMPL*)result; //1024
    fp1024_nayuki.execute_reverse_int(r->coefsC, p);
}
EXPORT void TorusPolynomial_ifft(LagrangeHalfCPolynomial* result, const TorusPolynomial* p) {
    LagrangeHalfCPolynomial_IMPL* r = (LagrangeHalfCPolynomial_IMPL*)result;
    fp1024_nayuki.execute_reverse_torus32(r->coefsC, p);
}
EXPORT void TorusPolynomial_fft(TorusPolynomial* result, const LagrangeHalfCPolynomial* p) {
    LagrangeHalfCPolynomial_IMPL* r = (LagrangeHalfCPolynomial_IMPL*)p;
    fp1024_nayuki.execute_direct_torus32(result, r->coefsC);
}







EXPORT uint64_t mod_mul(uint64_t x, uint64_t y, uint64_t p) {


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




}

//////////////////////////////////////////////////////////////////////////

EXPORT void bitrev_shuffle(uint64_t x[], int N) {

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

EXPORT void ntt_CT(uint64_t x[], uint64_t modulus, int N) {
    uint64_t b = 0;

    bitrev_shuffle(x, N);

    int wbarr = 0;

    for (int trans_size = 2; trans_size <= N; trans_size = trans_size * 2) {
        uint64_t wb = 1;
        for (int t = 0; t < (trans_size >> 1); t++) {
            for (int trans = 0; trans < (N / trans_size); trans++) {
                int i = trans * trans_size + t;
                int j = i + (trans_size >> 1);
                uint64_t a = x[i];
                if (t == 0) {
                    b = x[j];
                }
                else {

                    b = mod_mul(x[j], wb, modulus);
                }
                x[i] = mod_add(a, b, modulus);
                x[j] = mod_sub(a, b, modulus);
            }
            wb = wb_normal[wbarr];
            wbarr++;
        }
    }
}

EXPORT void intt_CT(uint64_t x[], uint64_t modulus, uint64_t y[], int N) {

    uint64_t b = 0;
    int inv = 0;

    bitrev_shuffle(x, N);

    for (int trans_size = 2; trans_size <= N; trans_size = trans_size * 2) {
        uint64_t wb = 1;

        for (int t = 0; t < (trans_size >> 1); t++) {
            for (int trans = 0; trans < (N / trans_size); trans++) {
                int i = trans * trans_size + t;
                int j = i + (trans_size >> 1);
                uint64_t a = x[i];
                if (t == 0)
                    b = x[j];
                else
                    b = mod_mul(x[j], wb, modulus);

                x[i] = mod_add(a, b, modulus);
                x[j] = mod_sub(a, b, modulus);
            }
            wb = wb_inverse[inv];
            inv++;
        }
    }

    int i = 0;
    y[0] = mod_mul(x[0], scale, modulus);

    for (i = 1; i < N; i++) {
        y[i] = mod_mul(x[i], scale, modulus); // scale*phi inversev, modulus  (phi inverse sclaed)
        uint64_t phi = phi_inverse[i];
        y[i] = mod_mul(y[i], phi, modulus);
    }
}