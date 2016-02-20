#ifndef CONTROLLER_CSC_H
#define CONTROLLER_CSC_H

#include "csntc.h"

#include <fstream>
#include <armadillo>
using namespace arma;

class Controller 
{

    public:

        Controller(int s, int n); 
        ~Controller(){}; 

        const vec &step(int t, int tc, int ti, double d, bool learn);

        const int trials = 3;

        int n_readout;

        shared_ptr<CSNTC> csntc;

        vec inp;
        vec inp_prev;
        vec activity;
        vec activity_prev;
        vec err;
        vec err_prev;

        vec readout_pot;
        vec readout;

        mat W_readout;
        vec W_env2CRX;

        ofstream trialstream;

    private:

        Controller(){};

};

#endif //CONTROLLER_CSC_H
