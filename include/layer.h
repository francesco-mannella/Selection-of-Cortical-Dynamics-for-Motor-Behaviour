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

#ifndef LAYER_H
#define LAYER_H

#include <parametrizable.h>
#include <armadillo>
#include <utils.h>

using namespace std;
using namespace arma;

/**
 * @brief Simple layer of leaky-integrator units
 */
class Layer : public Parametrizable
{
 
    public:
    
        Layer(std::string name); 
       
        Layer() = delete;
        Layer(const Layer &) = delete;
        Layer &operator=(const Layer &) = delete;

        virtual ~Layer(){};
       
        ////////////////////////////////////

        /** 
         * @brief Spreading step
         *
         * @param inp  input vector
         */
        virtual void step(vec inp);  
        
        /** 
         * @brief Store timestep activities into a timeline 
         *
         * @param t Timestep
         */     
        virtual void store(int t);
       
        /** 
         * @brief Reset activations 
         */
        virtual void reset();
        
        ////////////////////////////////////
       
        // VARS
        
        vec m_output; /*!< activation vector */
        mat m_data; /*!< activation storage */
        
        // CONSTS

        int m_STIME; /*!< simulation timesteps */
        double m_DT; /*!< integration step */
        double m_TAU; /*!< leaky dacay */
        double m_TH; /*!< leaky output threshold */
        double m_AMP; /*!< leaky output amplitude */
        int m_N; /*!< number of channels */
        
        double m_BL; /*!< activation baseline */
   
    private:
        
        // CONSTS
        
        
        vec m_u;
        vec m_v;
        
};


#endif //LAYER_H
