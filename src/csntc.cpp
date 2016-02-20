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

#include "csntc.h"

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


CSNTC::CSNTC(std::string name, int seed) :       
    Parametrizable(name),
    m_ganglia(name+"_ganglia"),
    m_thalamus(name+"_thalamus"),
    m_cortex(name+"_cortex", seed)
{
    
    m_cortex.init();
     
    addParameter("STIME",m_STIME);
    addParameter("THA2CRX_SPREAD",m_THA2CRX_SPREAD);
    addParameter("THA2CRX",m_THA2CRX);
    addParameter("CRX2THA",m_CRX2THA);
    addParameter("CRX2STR",m_CRX2STR);
    addParameter("GPI2THA",m_GPI2THA);

    try
    {
        loadParameters();
    }
    catch(parameter_file_exception &e)
    {
       
        m_STIME = 4000;

        m_THA2CRX_SPREAD = 0.1;
        m_THA2CRX = 1.5;
        m_CRX2THA = 0.3;
        m_CRX2STR = 0.1;
        m_GPI2THA = 20.0;

        saveParameters();
    }

    m_N_CHA = m_ganglia.m_N;
    m_N_THA = m_thalamus.m_N;
    m_N_CRX = m_cortex.m_N;


    // number of CRX cells connected to each GPI channel
    int T2CP = ceil(m_N_CRX*m_THA2CRX_SPREAD);  
    
    // relative indices of connected CRX cells ()
    uvec rng = conv_to<uvec>::from(linspace(0,T2CP-1,T2CP));
    
    // unary zero vector (for armadillo range operations)
    uvec elem = conv_to<uvec>::from(zeros(1));

    m_GPI2THA_w = zeros<vec>(m_N_CHA);
    m_THA2CRX_W = zeros<mat>(m_N_CRX,m_N_THA);
    m_CRX2THA_W = zeros<mat>(m_N_THA,m_N_CRX);
    m_CRX2STR_W = zeros<mat>(m_N_CHA,m_N_CRX);

    // weights from GPI to THA are negative
    {
        vec &w = m_GPI2THA_w;

        w = -ones(m_N_CHA)*m_GPI2THA;
    }

    // each THA channel reaches only the relative CRX 
    // subpopulation
    {
        mat &W = m_THA2CRX_W;

        for( int x=0; x<m_N_THA; x++) 
            W( rng + T2CP*x , elem+x ).fill( m_THA2CRX);
    }

    // units in each CRX subpopulation reach only the relative
    // THA channel
    {
        mat &W = m_CRX2THA_W;

        for( int x=0; x<m_N_THA; x++) 
            W( elem+x , rng + T2CP*x ).fill(m_CRX2THA);
    }

    // units in each CRX subpopulation reach only the relative
    // STR channel
    {
        mat &W = m_CRX2STR_W;

        for( int x=0; x<m_N_CHA; x++) 
            W( elem+x , rng + T2CP*x ).fill(m_CRX2STR/T2CP);
    }

    // initialize all variables
    reset();

    // initialize the data storage array
    m_data = zeros<mat>( 
            m_ganglia.m_N + 
            m_thalamus.m_N + 
            m_cortex.m_N, m_STIME );
}

////////////////////////////////////////////////////////////////////////////

void CSNTC::reset()
{
    m_ganglia.reset();
    m_thalamus.reset();
    m_cortex.reset();
}

////////////////////////////////////////////////////////////////////////////

void CSNTC::step(vec crx, vec str, vec da)
{
        // ganglia input
        vec bg_str_inp = str;
        vec bg_crx_inp = m_CRX2STR_W*m_cortex.m_output;
        
        // thalamus input
        vec th_inp =
            m_GPI2THA_w%m_ganglia.m_gpi + 
            m_CRX2THA_W*m_cortex.m_output;

        // cortex input
        vec crx_inp = 
            crx + 
            m_THA2CRX_W*m_thalamus.m_output;


        
        m_ganglia.step(bg_str_inp, bg_crx_inp, da);
        m_thalamus.step(th_inp);
        m_cortex.step(crx_inp);


}


////////////////////////////////////////////////////////////////////////////

void CSNTC::store(int t)
{
    m_ganglia.store(t);
    m_thalamus.store(t);
    m_cortex.store(t);
}

////////////////////////////////////////////////////////////////////////////

