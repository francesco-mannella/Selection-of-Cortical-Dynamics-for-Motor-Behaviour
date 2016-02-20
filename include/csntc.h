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

#ifndef CSNTC_H
#define CSNTC_H

/**
 * @mainpage CSNTC - A model of the Cortico-Striato-Nigro-Thalamo-Cortical loop
 *
 * This library gives an implementation of the CSNTC loops by putting together 
 * a biologically plausible model of hte basal ganglia (BasalGanglia),
 * An echo-state network as a model of cortex (Reservoir) and a simplified 
 * representation of the thalamus (Layer).
 *
 *
 */


#include <parametrizable.h>
#include <layer.h>
#include <basalganglia.h>
#include <reservoir.h>
#include <utils.h>


#include <armadillo>

using namespace std;
using namespace arma;

/**
 * @brief Cortico-striato-nigro-thalamo-cortical module
 */
class CSNTC : public Parametrizable
{
 
    public:
    
        CSNTC(std::string name, int seed = 1); 
        
        CSNTC() = delete;
        CSNTC(const CSNTC &) = delete;
        CSNTC &operator=(const CSNTC &) = delete;
        
        virtual ~CSNTC(){};
       
        /** 
         * @brief Spreading step
         *
         * @param crx   Reentrant cortical input vector
         * @param str   Striatal External input vector
         * @param da     Dopamine 
         */
        virtual void step( vec crx, vec str, vec da); 
 
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
         
        BasalGanglia m_ganglia;    /*!< Basal ganglia module */
        Layer m_thalamus;    /*!< Thalamic module */
        Reservoir m_cortex;    /*!< Cortical module */
        
        mat m_data;    /*!< Store of activation hystory */
        
        // CONSTS

        int m_STIME;    /*!< Simulation timesteps */
        
        double m_N_CRX;    /*!< Number of cortical units */
        double m_N_CHA;    /*!< Number of basal ganglia channels */
        double m_N_THA;    /*!< Number of thalamic units */ 
        
        // CONSTS

        double m_THA2CRX_SPREAD;    /*!< Spread of thalamic input 
                                         over the cortical units */
        double m_THA2CRX;    /*!< Thalamus -> Cortex weight parameter */
        double m_CRX2THA;    /*!< Cortex -> Thalamus weight parameter */
        double m_CRX2STR;    /*!< Cortex -> Striatum weight parameter */
        double m_GPI2THA;    /*!< GPi -> Thalamus weight parameter */
        
        mat m_THA2CRX_W;    /*!< Thalamus -> Cortex */
        mat m_CRX2THA_W;    /*!< Cortex -> Thalamus */
        mat m_CRX2STR_W;    /*!< Cortex -> Striatum */
        vec m_GPI2THA_w;    /*!< GPi -> Thalamus */


};

#endif // CSNTC_H
