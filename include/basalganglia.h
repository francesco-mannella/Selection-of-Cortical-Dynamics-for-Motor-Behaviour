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

#ifndef BASAL_GANGLIA_H
#define BASAL_GANGLIA_H

#include <parametrizable.h>
#include <utils.h>
#include <armadillo>

using namespace std;
using namespace arma;

/**
 * @brief Basal ganglia module   
 *
 * This class implements the computational hypothesis that
 * the basal ganglia act as a selector due to their 
 * intrinsic competitive mechanisms.
 *
 *
 */
class BasalGanglia : public Parametrizable
{
 
    public:
    
        
        BasalGanglia(std::string name); 
        
        BasalGanglia() = delete;
        BasalGanglia(const BasalGanglia &) = delete;
        BasalGanglia &operator=(const BasalGanglia &) = delete;

        virtual ~BasalGanglia(){};
       
        ////////////////////////////////////
        
        /** 
         * @brief Spreading step
         *
         * @param inp External input vector
         * @param crx Reentrant cortical input vector
         * @param da  Dopamine 
         */
        virtual void step( vec inp, vec crx, vec da);      
         
        /** 
         * @brief Store timestep activities into a timeline 
         *
         * @param t Timestep
         */       
        virtual void store(int t);
        
        /**
         * @brief Save timeline on file 
         */
        virtual void save();    
        
        /** 
         * @brief Reset activations 
         */
        virtual void reset();          

        ////////////////////////////////////
       
        // VARS

        // Activations

        vec m_sd1;    /*!< Striatum D1 - activation */ 
        vec m_sd2;    /*!< Striatum D2 - activation  */
        vec m_stn;    /*!< Subthalamic nucleus - activation  */
        vec m_gpi;    /*!< Internal globus pallidus - activation  */
        vec m_gpe;    /*!< External globus pallidus - activation  */
 
        // Potentials

        vec m_sd1p;    /*!< Striatum D1 - potential */ 
        vec m_sd2p;    /*!< Striatum D2 - potential  */
        vec m_stnp;    /*!< Subthalamic nucleus - potential  */
        vec m_gpip;    /*!< Internal globus pallidus - potential  */
        vec m_gpep;    /*!< External globus pallidus - potential  */

        // Storage

        mat m_data;    /*!< Storage of the timeline of activations */
        
        // CONSTS

        int m_STIME;    /*!< Simulation timesteps */
        double m_DT;    /*!< Integration step */
        double m_TAU;    /*!< Leaky decay */
        double m_TH;    /*!< Leaky output threshold */
        double m_AMP;    /*!< Leaky output amplitude */
        int m_N;    /*!< Number of channels */
        static const int m_NUMCOMP = 5;     /*!< Number of layers */ 
 
        
        
        double m_SD1_BL;    /*!< Baseline activation of striatum D1 */
        double m_SD1_DA;    /*!< DA multiplying factor of striatum D1 */
        double m_SD2_BL;    /*!< Baseline activation of striatum D2 */
        double m_SD2_DA;    /*!< DA multiplying factor of striatum D2 */
        double m_GPI_BL;    /*!< Baseline activation of GPi */
        double m_GPE_BL;    /*!< Baseline activation of striatum of GPe */

        
        // WEIGHTS

        double m_INP2SD1_W;    /*!< Input -> striatum D1 */ 
        double m_INP2SD2_W;    /*!< Input -> striatum D2 */
        double m_INP2STN_W;    /*!< Input -> STN */
        double m_CRX2SD1_W;    /*!< Cortex -> striatum D1 */
        double m_CRX2SD2_W;    /*!< Cortex -> striatum D2*/
        double m_CRX2STN_W;    /*!< Cortex -> STN */
        double m_SD12GPI_W;    /*!< Striatum D1 -> GPi */
        double m_SD22GPE_W;    /*!< Striatum D2 -> GPe */
        double m_STN2GPI_W;    /*!< STN -> GPi */
        double m_STN2GPE_W;    /*!< STN -> GPe */
        double m_GPE2STN_W;    /*!< GPe -> STN */
        double m_GPE2GPI_W;    /*!< GPe -> GPi */
        
};

#endif //BASAL_GANGLIA_H
