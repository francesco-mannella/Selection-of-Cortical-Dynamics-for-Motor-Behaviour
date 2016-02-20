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

#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <csntc.h>
#include <iostream>
#include <utils.h>
#include <memory>


/**
 * @brief Manager of simulations in Figs 7 and 8
 *
 * A class containing all data and methods to run a simulation. 
 * Inherits methods to synchronize parameters with a configuration file from Parametrizable. 
 */
struct Simulator : public Parametrizable
{

    //////////////////////////////////////////////////////////////////////////////////////////////
    //// Methods /////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @param   name        A string defining the name of the simulator object (used for parameter storage)
     * @param   SEED        The seed for random number generation 
     */
    Simulator(string name, int SEED = 1);
   
    // Cannot use empty, copy or assign constructors
    Simulator() = delete;
    Simulator(const Simulator& copy) = delete;
    Simulator& operator=(const Simulator& copy) = delete;

    ~Simulator(){};

    void build_funcs(const string &teach_sec_file, 
            const string &fun_sec_file = "", bool reset_all=false);     ///< Build input and theach time series 
    void modify_trial_sequence(bool reset=false);    ///< Shuffle per-trial timeseries of input an teaching of the single session 
    void run_offline();    ///< Run a life-cicle for offline learning simulation 
    void run_offline_double_teach();    ///< Run a life-cicle for offline learning simulation 
    void run_online();    ///< Run a life-cicle for offline learning simulation 
    void step(int t);    ///< Spreading step 
    void step_noise(int t);    ///< Spreading step with further noise in the striatal input (for striato-cortical learning) 

    //////////////////////////////////////////////////////////////////////////////////////////////
    //// Members /////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////
    
    // Members initialized from file

    int m_SESSIONS;    ///< Number of sessions - ONLINE
    int m_STIME;    ///< Number of total timesteps per session
    int m_TTIME;    ///< Number of timesteps per trial
    int m_SELECTIONS;    ///< Number of possible selections (bg channels)

    double m_FUN_PERIOD;    ///< Period of the cortical input wave function 
    double m_FUN_CENTER;    ///< Equilibrium center of the cortical input wave function
    double m_FUN_PHASE;    ///< Phase of the cortical input wave function
    double m_FUN_AMP;    ///< Amplitude of the cortical input wave function

    double m_BG_BIAS_AMP;    ///< Bias amplitude of the input to the to-be-selected striatal channel
    double m_BG_NOISE;    ///< Noise to the striatum 

    double m_DT;    ///< Integration timestep
    double m_TAU;    ///< Leaky decay
    double m_DELTA;    ///< DT/TAU
    double m_LAMBDA;    ///< Ridge regression regularization parameter - OFFLINE
    double m_BETA;    ///< BPDC regularization parameter - ONLINE
    double m_ETA;    ///< Learning rate - ONLINE
   
    double m_LEARN_WINDOW;    ///< Length of the learning window in timesteps

    double m_DA_GAP;    ///< Time gap of dopaminergic switch-off at start and end of trial
    double m_INP2M;     ///< Maximum amplitude of weights from input to motor cortex
    double m_INP2BG;     ///< Maximum amplitude of weights from input to striatum

    string m_PRINT;    ///< "TRUE" full storage of data of the model - "FALSE" otherwise


    // Simulation members 
    unique_ptr<CSNTC> m_csntc;    ///< A CSNTC object containing all data and methods for a CSNTC module
    // Readout
    mat m_rout;    ///< Read-out story 
    mat m_Wout;    ///< Read-out weights 
    int m_n_readout;    ///< Number of readouts
    // Weights
    double m_W_inp2bg;    ///< Read-out  weights
    mat m_W_inp2m;    ///< Read-out  weights
    // Inputs
    vec m_da;    ///< A vector of length m_STIME containing the dopaminergic modulation    
    mat m_fun;    ///< The cortical wave input during a session (length = m_STIME) 
    mat m_crx_inp;    ///< Vector of input to cortical units 
    mat m_bg_inp; ///< Vector of input to the striatum

    // RANGES
    //uvec m_rng1;   ///< Timestep m_crx_idx of the learning window within the 1-st trial of a session (index 0 is start of session)
    //uvec m_rng2;   ///< Timestep m_crx_idx of the learning window within the 2-nd trial of a session
    //uvec m_rng3;   ///< Timestep m_crx_idx of the learning window within the 3-rd trial of a session
    uvec m_rng;    ///< [ m_rng1' m_rng2' m_rng3' ] - matlab syntax
    umat m_rngs;    ///< [ m_rng1 m_rng2 m_rng3 ]' - matlab syntax  

    // TEACH
    mat m_teach;    ///< Each row contains a s_STIME vector with the target function of a readout unit for each trial 
    mat m_curr_teach;    ///< Same as m_teach except that the sequence of trials is shuffled on every session

    // Indexes for permutations
    uvec m_crx_idx;    ///< Indices of the vector of cortical units
    uvec m_trial_idx;   ///< Trial sequence

    static const string m_TEACH_FILE_1;    ///< Filename for data about target activities
    static const string m_TEACH_FILE_2;    ///< Filename for data about target activities
    static const string m_TEACH_FILE_3;    ///< Filename for data about target activities
    static const string m_TEACH_FILE_4;    ///< Filename for data about target activities
    static const string m_FUN_FILE_1;    ///< Filename for data about input 
    static const string m_FUN_FILE_2;    ///< Filename for data about input 
    static const string m_FUN_FILE_3;    ///< Filename for data about input 
    static const string m_FUN_FILE_4;    ///< Filename for data about input 

};


#endif //SIMULATOR_H

