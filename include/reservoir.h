/* * CSNTC - Cortico-striato-nigro-thalamo-cortical comutational model
 * Copyright (C) 2014 Francesco Mannella <francesco.mannella@gmail.com>
 *
 * This file is part of CSNTC.
 *
 * CSNTC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CSNTC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CSNTC.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef RESERVOIR_H
#define RESERVOIR_H

#include <iostream>
#include <armadillo>
#include <parametrizable.h>
#include <utils.h>

using namespace std;
using namespace arma;


/**
 * @brief Firing-rate echo-state network
 */
class Reservoir : public Parametrizable 
{

    public :

        // METHODS

        Reservoir(string file,int seed=1);
          
        Reservoir() = delete;
        Reservoir(const Reservoir &) = delete;
        Reservoir &operator=(const Reservoir &) = delete;

       ~Reservoir(){}
   
        /** initialize vectors */
        void init();     

        /**
         * normalize weights to spectral radius is rho
         *
         * @param w    a weight matrix
         * @param rho  the desired spectral radius
         * @return     the modified weight matrix
         */ 
        mat normalize(mat w, double rho);   

        /**
         * make a random sparse matrix with modulated 
         * infinitesimal rotation and translation 
         *
         * @return a weight matrix
         */ 
        mat randomize();    

        /** 
         * find the optimized spectral radius of W so that 
         * rho(1-epsilon) < Wd  < 1, where Wd = (dt/tau)*W+(1-(dt/tau))*eye. 
         * See Proposition 2 in Jaeger et al. (2007) http://goo.gl/bqGAJu.
          */
        void normalize_to_echo();     
        
        
        void reset();    /*!< reset activations */
        void store(int t);    /*!< store */
        void step(vec inp);    /*!< single step */

        // CONSTS

        string m_NAME;    /*!< name of the object */
        double m_STIME;    /*!< Simulation timesteps */
        double m_DT;    /*!< Integration step */
        double m_TAU;    /*!< Leaky decay */
        double m_TH;    /*!< Leaky output threshold */
        double m_AMP;    /*!< Leaky output amplitude */
        double m_N;    /*!< Number of channels */
        bool m_TRUNK;   /*!< if tanh output function is trunked */
        bool m_NOISE;    /*!< switch on/off internal noise */
        double m_NOISE_STD;    /*!<internal noise standard deviation */
        string m_sTRUNK;    /*!< string version of m_TRUNK for file sync */
        double m_SPARSENESS;    /*!< sparseness proportion of internal connections */
        double m_WEIGHT_MEAN;    /*!< mean value of weights */
        string m_sNOISE;    /*|< string version of m_NOISE for file sync*/ 
        double m_ALPHA;    /*!< traslation eigenvalues proportion */
        double m_BETA;    /*!< rotation eigenvalues proportion */
        double m_GAMMA;    /*!<  mean point of eigenvalues */
        double m_EPSILON;    /*!<  distance from the spectral radius limit for chaotic behaviour */
        double m_G; /*!<  spectral amplification */

        int m_SEED;   /*!< seed of the random generation (for file sync) */
        string m_LOAD;    /*!< switch on/off loading weights from file */
 
        // VARS

        vec m_output;    /*!< activaion output */
        mat m_w;    /*! internal weights */
        mat m_data;    /*!< activaion storage */

    protected:

        vec m_u;
        vec m_v;

};

#endif //RESERVOIR_H

