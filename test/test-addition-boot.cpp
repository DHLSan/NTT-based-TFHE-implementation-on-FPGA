#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "tfhe.h"
#include "polynomials.h"
#include "lwesamples.h"
#include "lwekey.h"
#include "lweparams.h"
#include "tlwe.h"
#include "tgsw.h"
#include <time.h>

using namespace std;



// **********************************************************************************
// ********************************* MAIN *******************************************
// **********************************************************************************


void dieDramatically(string message) {
    cerr << message << endl;
    abort();
}


#ifndef NDEBUG
extern const TLweKey* debug_accum_key;
extern const LweKey* debug_extract_key;
extern const LweKey* debug_in_key;
#endif

int32_t main(int32_t argc, char** argv) {
#ifndef NDEBUG
    cout << "DEBUG MODE!" << endl;
#endif
    const int32_t nb_bits = 64;
    const int32_t nb_trials = 1;


    int arrm1[] = { 1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }; //3
    int arrm2[] = { 0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }; //26


    // generate params
    int32_t minimum_lambda = 70;
    TFheGateBootstrappingParameterSet* params = new_default_gate_bootstrapping_parameters(minimum_lambda);
    const LweParams* in_out_params = params->in_out_params;

    // generate the secret keyset
    TFheGateBootstrappingSecretKeySet* keyset = new_random_gate_bootstrapping_secret_keyset(params);

    // Give input outputs
    cout << "First message : ";
    for (int i = 0; i < nb_bits; i++) {
        cout << arrm1[i];
    }
    cout << endl;
    cout << "Second message : ";
    for (int i = 0; i < nb_bits; i++) {
        cout << arrm2[i];
    }

    cout << "" << endl;
    cout << "" << endl;



    for (int32_t trial = 0; trial < nb_trials; ++trial) {


        LweSample* x = new_LweSample_array(nb_bits, in_out_params);
        LweSample* y = new_LweSample_array(nb_bits, in_out_params);
        LweSample* temp = new_LweSample_array(nb_bits, in_out_params);

        for (int32_t i = 0; i < nb_bits; ++i) {
            bootsSymEncrypt(x + i, arrm1[i], keyset);
            bootsSymEncrypt(y + i, arrm2[i], keyset);
        }  // ONE AND ZERO ARRAY encrypting



        //*************************************************XOR GATE**********************************************************
        cout << "The result of Homomorphic OR:   ";
        for (int i = 0; i < nb_bits; i++) {
            bootsOR(temp + i, x + i, y + i, &keyset->cloud);
            bool messSum = bootsSymDecrypt(temp + i, keyset);
            cout << messSum;
        }

   
      

    }

    delete_gate_bootstrapping_secret_keyset(keyset);
    delete_gate_bootstrapping_parameters(params);

    return 0;
}