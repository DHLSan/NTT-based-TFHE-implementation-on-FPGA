#include <polynomials.h>
#include "lagrangehalfc_impl.h"
#include <iostream>
using namespace std;
#define mod_add(a, b, modulus) ((modulus - a) > b) ? (a + b) : (a + b - modulus)
#define mod_sub(a, b, modulus) (a >= b) ? (a - b) : (modulus - b + a)
LagrangeHalfCPolynomial_IMPL::LagrangeHalfCPolynomial_IMPL(const int32_t N) {

    coefsC = new uint64_t[1024];
    proc = &fp1024_nayuki;
}

LagrangeHalfCPolynomial_IMPL::~LagrangeHalfCPolynomial_IMPL() {
    delete[] coefsC;
}

//initialize the key structure
//(equivalent of the C++ constructor)
EXPORT void init_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* obj, const int32_t N) {
    new(obj) LagrangeHalfCPolynomial_IMPL(N);
}
EXPORT void init_LagrangeHalfCPolynomial_array(int32_t nbelts, LagrangeHalfCPolynomial* obj, const int32_t N) {
    for (int32_t i = 0; i < nbelts; i++) {
        new(obj + i) LagrangeHalfCPolynomial_IMPL(N);
    }
}

//destroys the LagrangeHalfCPolynomial structure
//(equivalent of the C++ destructor)
EXPORT void destroy_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* obj) {
    LagrangeHalfCPolynomial_IMPL* objbis = (LagrangeHalfCPolynomial_IMPL*)obj;
    objbis->~LagrangeHalfCPolynomial_IMPL();
}
EXPORT void destroy_LagrangeHalfCPolynomial_array(int32_t nbelts, LagrangeHalfCPolynomial* obj) {
    LagrangeHalfCPolynomial_IMPL* objbis = (LagrangeHalfCPolynomial_IMPL*)obj;
    for (int32_t i = 0; i < nbelts; i++) {
        (objbis + i)->~LagrangeHalfCPolynomial_IMPL();
    }
}


//MISC OPERATIONS
/** sets to zero */
EXPORT void LagrangeHalfCPolynomialClear(
    LagrangeHalfCPolynomial* reps) {
    LagrangeHalfCPolynomial_IMPL* reps1 = (LagrangeHalfCPolynomial_IMPL*)reps;
    const int32_t Ns2 = 1024;
    for (int32_t i = 0; i < 1024; i++)
        reps1->coefsC[i] = 0;
}

EXPORT void LagrangeHalfCPolynomialSetTorusConstant(LagrangeHalfCPolynomial* result, const Torus32 mu) {
    LagrangeHalfCPolynomial_IMPL* result1 = (LagrangeHalfCPolynomial_IMPL*)result;
    const int32_t Ns2 = result1->proc->Ns2;
    uint64_t* b = result1->coefsC;

    for (int32_t j = 0; j < Ns2; j++)
        b[j] = (uint64_t)mu;


}

EXPORT void LagrangeHalfCPolynomialAddTorusConstant(LagrangeHalfCPolynomial* result, const Torus32 mu) {
    LagrangeHalfCPolynomial_IMPL* result1 = (LagrangeHalfCPolynomial_IMPL*)result;
    const int32_t Ns2 = result1->proc->Ns2;
    uint64_t* b = result1->coefsC;

    for (int32_t j = 0; j < Ns2; j++)
        b[j] += mu;

}

EXPORT void LagrangeHalfCPolynomialSetXaiMinusOne(LagrangeHalfCPolynomial* result, const int32_t ai) {

}

/** termwise multiplication in Lagrange space */
EXPORT void LagrangeHalfCPolynomialMul(
    LagrangeHalfCPolynomial* result,
    const LagrangeHalfCPolynomial* a,
    const LagrangeHalfCPolynomial* b) {
    LagrangeHalfCPolynomial_IMPL* result1 = (LagrangeHalfCPolynomial_IMPL*)result;
    const int32_t Ns2 = result1->proc->Ns2;
    uint64_t* aa = ((LagrangeHalfCPolynomial_IMPL*)a)->coefsC;
    uint64_t* bb = ((LagrangeHalfCPolynomial_IMPL*)b)->coefsC;
    uint64_t* rr = result1->coefsC;
    for (int32_t i = 0; i < 1024; i++)
        rr[i] = mod_mul(aa[i], bb[i], 0xffffffff00000001);
}

/** termwise multiplication and addTo in Lagrange space */
EXPORT void LagrangeHalfCPolynomialAddMul(
    LagrangeHalfCPolynomial* accum,
    const LagrangeHalfCPolynomial* a,
    const LagrangeHalfCPolynomial* b)
{
    LagrangeHalfCPolynomial_IMPL* result1 = (LagrangeHalfCPolynomial_IMPL*)accum;
    const int32_t Ns2 = result1->proc->Ns2;
    uint64_t* aa = ((LagrangeHalfCPolynomial_IMPL*)a)->coefsC;
    uint64_t* bb = ((LagrangeHalfCPolynomial_IMPL*)b)->coefsC;
    uint64_t* rr = result1->coefsC;
    uint64_t temp[1024];
    for (int32_t i = 0; i < 1024; i++) {
        temp[i] = mod_mul(aa[i], bb[i], 0xffffffff00000001); // temp
        rr[i] = mod_add(rr[i], temp[i], 0xffffffff00000001);
    }
}


/** termwise multiplication and addTo in Lagrange space */
EXPORT void LagrangeHalfCPolynomialSubMul(
    LagrangeHalfCPolynomial* accum,
    const LagrangeHalfCPolynomial* a,
    const LagrangeHalfCPolynomial* b)
{
    LagrangeHalfCPolynomial_IMPL* result1 = (LagrangeHalfCPolynomial_IMPL*)accum;
    const int32_t Ns2 = result1->proc->Ns2;
    uint64_t* aa = ((LagrangeHalfCPolynomial_IMPL*)a)->coefsC;
    uint64_t* bb = ((LagrangeHalfCPolynomial_IMPL*)b)->coefsC;
    uint64_t* rr = result1->coefsC;
    uint64_t temp[1024];
    for (int32_t i = 0; i < 1024; i++) {
        temp[i] = mod_mul(aa[i], bb[i], 0xffffffff00000001);
        rr[i] = mod_sub(rr[i], temp[i], 0xffffffff00000001);
    }
}

EXPORT void LagrangeHalfCPolynomialAddTo(
    LagrangeHalfCPolynomial* accum,
    const LagrangeHalfCPolynomial* a) {
    LagrangeHalfCPolynomial_IMPL* result1 = (LagrangeHalfCPolynomial_IMPL*)accum;
    const int32_t Ns2 = result1->proc->Ns2;
    uint64_t* aa = ((LagrangeHalfCPolynomial_IMPL*)a)->coefsC;
    uint64_t* rr = result1->coefsC;
    for (int32_t i = 0; i < 1024; i++) {
        rr[i] = mod_add(aa[i], rr[i], 0xffffffff00000001);
    }
}

