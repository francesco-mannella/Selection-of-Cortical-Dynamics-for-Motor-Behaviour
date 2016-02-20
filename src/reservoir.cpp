//
// CSNTC - Cortico-striato-nigro-thalamo-cortical comutational model
// Copyright (C) 2014 Francesco Mannella <francesco.mannella@gmail.com>
//
// This file is part of CSNTC.
//
// CSNTC is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// CSNTC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with CSNTC.  If not, see <http://www.gnu.org/licenses/>.
//

#include <reservoir.h>
#include <md5.h>
#include <iomanip>
#include <stdexcept>

/**
 * Evaluate the spectral radius of a square matrix
 *
 * @param    m     a matrix
 * @return   the spectral radius
 *   
 */
double evaluate_rho(mat m)
{
    if (m.n_cols != m.n_rows)
        throw std::invalid_argument( "evaluate_rho(m) requires square matrix" );
    
    // find eigenvalues
    cx_vec eigenvalues; 
    cx_mat eigenvectors; 
    eig_gen(eigenvalues, eigenvectors, m);
    
    // the spectral radius is the maximum between the 
    // absolute values of eigenvalues.
    return  abs(eigenvalues).max();

}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

Reservoir::Reservoir(string file,int seed) : Parametrizable(file) 
{
    m_SEED = seed;
}

////////////////////////////////////////////////////////////////////////////

void Reservoir::init()
{

    addParameter("DT",m_DT);
    addParameter("TAU",m_TAU);
    addParameter("STIME",m_STIME);
    addParameter("N",m_N);
    addParameter("ALPHA",m_ALPHA);
    addParameter("BETA",m_BETA);
    addParameter("GAMMA",m_GAMMA);
    addParameter("TH",m_TH);
    addParameter("AMP",m_AMP);
    addParameter("TRUNK",m_sTRUNK);
    addParameter("NOISE",m_sNOISE);
    addParameter("NOISE_STD",m_NOISE_STD);
    addParameter("G",m_G);
    addParameter("EPSILON",m_EPSILON);
    addParameter("SPARSENESS",m_SPARSENESS);
    addParameter("WEIGHT_MEAN",m_WEIGHT_MEAN);
    addParameter("LOAD_WEIGHTS",m_LOAD);

    try
    {
        loadParameters();
    }
    catch(parameter_file_exception &e)
    {
       
        m_DT = 0.001;
        m_TAU = 0.004;
        m_STIME = 400;
        m_N = 150;
        m_ALPHA = 0.5;
        m_BETA = 0.5;
        m_GAMMA = 0.0;
        m_TH = 0.0;
        m_AMP = 1.0;
        m_sNOISE = "FALSE";
        m_NOISE_STD = 0.1;
        m_G = 1.0;
        m_EPSILON = 0.0001;
        m_SPARSENESS = 1.0;
        m_WEIGHT_MEAN = 0.0;
        m_sTRUNK = "FALSE";
        m_LOAD = "FALSE";
        
        saveParameters();

    }


    m_TRUNK = (m_sTRUNK=="TRUE");
    m_NOISE = (m_sNOISE=="TRUE");

    m_u = zeros<vec>(m_N);
    m_v = m_u;
    m_output = m_u;
          
   
    // use md5 hash to get the name of the weights file
    stringstream ss;
    ss << m_filename << m_SEED << 
        m_DT << m_TAU << m_STIME << 
        m_N << m_ALPHA<< m_BETA << 
        m_GAMMA << m_TH<< m_AMP<< 
        m_sNOISE << m_NOISE_STD<< m_G << 
        m_EPSILON<< m_SPARSENESS << m_WEIGHT_MEAN; 
    string fweights = m_filename  + "_" + md5(ss.str());

    // load/save weights from/to file
    if(m_LOAD=="TRUE")
    {
        if(ifstream(fweights.c_str()).good() )
        {
            m_w.load(fweights);
        }
        else
        {
            normalize_to_echo();
            m_w.save(fweights);
            m_w.load(fweights);
        }
    }
    // initialize weights
    else
    {
        normalize_to_echo();
    }

    // initialize data storage array
    m_data = zeros<mat>(m_N,m_STIME);
}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

void Reservoir::reset()
{
    m_u *= 0;
    m_v *= 0;
    m_output *= 0;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

void Reservoir::step(vec inp)
{
    // unit integration
    m_u += (m_DT/m_TAU)*( -m_u + inp + m_w*m_v);
    if (m_NOISE ==true)
        m_u += (m_DT/m_TAU)*(m_NOISE_STD*randn(m_N) );
    m_v = outfun( m_u, m_TH, m_AMP);
       
    // truncated activation function
    if (m_TRUNK == true) 
        m_v = m_v%(m_v>0);
    
    m_output = m_v;
}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

void Reservoir::store(int t)
{
    m_data(span::all,t) = m_output;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

mat Reservoir::randomize()
{
    // initial random choise + sparseness
    mat w = (randn<mat>(m_N,m_N) + m_WEIGHT_MEAN)%
        (randu<mat>(m_N,m_N)<m_SPARSENESS);
   
    // remodulate rotation expansion,
    // contraction/espansion and eigenvalues average
    w = 
        m_ALPHA * (w + w.t()) + 
        m_BETA * (w - w.t()) +
        m_GAMMA * eye(m_N,m_N);

    return w;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

mat Reservoir::normalize(mat w, double rho)
{ 
    // find spectral radius 
    double curr_rho = evaluate_rho(w);
    
    // move the spectral radius to rho
    return rho * w / curr_rho;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

void Reservoir::normalize_to_echo() 
{
    double target = 1-m_EPSILON/2.;    // we will try to reach this rho in the effective matrix
    double rho_estimation = 100.0;       // initial rho estimation
    double h = m_DT/m_TAU;   // convenience decay constant
    mat I = eye(m_N,m_N);    // convenience identity matrix
    
    // first create a matrix of rho 1 based on 
    // alpha beta and gammma parameters
    m_w = randomize();
    m_w = normalize(m_w, 1.0);

    // evaluate the effective dynamical matrix
    mat e_w = h*rho_estimation*m_w + (1 - h)*I;
  
    // initial effective rho estimation
    double effective_rho_estimation = evaluate_rho(e_w); 

    // we try to reach the leaky echo-state condition
    while (   not (( (1-m_EPSILON) < effective_rho_estimation) and (effective_rho_estimation<1.0)) ) 
    {
        // evaluate the effective dynamical matrix
        e_w = h*rho_estimation*m_w + (1 - h)*I;
     
        // evaluate its spectral radius
        effective_rho_estimation = evaluate_rho(e_w);
      
        // update the real spectral radius based on 
        // the amount of distance to the leaky-echo-state condition
        rho_estimation += (1.0/h)*(target-effective_rho_estimation);
    }
 
    // finally modify the weight matrix to the 
    // converged spectral radius
    m_w = rho_estimation*m_G*m_w;
 
}
