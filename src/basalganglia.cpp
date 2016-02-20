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

#include "basalganglia.h"

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


BasalGanglia::BasalGanglia(std::string name) :       
    Parametrizable(name)
{

    // Add members to the list of 
    // parameters to sync with file

    addParameter("STIME",m_STIME);
    addParameter("DT",m_DT);
    addParameter("TAU",m_TAU);
    addParameter("TH",m_TH);
    addParameter("AMP",m_AMP);
    addParameter("N",m_N);
    addParameter("INP2SD1_W",m_INP2SD1_W);
    addParameter("INP2SD2_W",m_INP2SD2_W);
    addParameter("INP2STN_W",m_INP2STN_W);
    addParameter("CRX2SD1_W",m_CRX2SD1_W);
    addParameter("CRX2SD2_W",m_CRX2SD2_W);
    addParameter("CRX2STN_W",m_CRX2STN_W);
    addParameter("SD12GPI_W",m_SD12GPI_W);
    addParameter("SD22GPE_W",m_SD22GPE_W);
    addParameter("STN2GPI_W",m_STN2GPI_W);
    addParameter("STN2GPE_W",m_STN2GPE_W);
    addParameter("GPE2STN_W",m_GPE2STN_W);
    addParameter("GPE2GPI_W",m_GPE2GPI_W); 
    addParameter("SD1_BL",m_SD1_BL);
    addParameter("SD1_DA",m_SD1_DA);
    addParameter("SD2_BL",m_SD2_BL);
    addParameter("SD2_DA",m_SD2_DA);
    addParameter("GPI_BL",m_GPI_BL);
    addParameter("GPE_BL",m_GPE_BL);

    // Try to load parameter file
    // if does not exist initialize 
    // parameters and create it for the 
    // next time

    try
    {
        loadParameters();
    }
    catch(parameter_file_exception &e)
    {

        m_STIME = 1000;
        m_DT = 0.001;
        m_TAU = 0.005;
        m_TH = 0.0;
        m_AMP = 1.0;
        m_N = 3;

        m_INP2SD1_W = 0.1;
        m_INP2SD2_W = 0.1;
        m_INP2STN_W = 0.0;
        m_CRX2SD1_W = 0.4;
        m_CRX2SD2_W = 3.2;
        m_CRX2STN_W = 1.0;
        m_SD12GPI_W = 3.0;
        m_SD22GPE_W = 2.0;
        m_STN2GPI_W = 1.0;
        m_STN2GPE_W = 1.0;
        m_GPE2STN_W = 0.8;
        m_GPE2GPI_W = 0.3;

        m_SD1_BL = 0.1;
        m_SD1_DA = 0.5;
        m_SD2_BL = 0.1;
        m_SD2_DA = 20.0;
        m_GPI_BL = 0.1;
        m_GPE_BL = 0.1;

        saveParameters();
    }

    // reset activations and storage
    
    reset();
    m_data = zeros<mat>( m_NUMCOMP*m_N, m_STIME );
}

////////////////////////////////////////////////////////////////////////////


void BasalGanglia::reset()
{
    
    // init all activation and potential 
    // vectors to zero
    
    m_sd1p = zeros<vec>(m_N);
    m_sd2p = zeros<vec>(m_N);
    m_stnp = zeros<vec>(m_N);
    m_gpip = zeros<vec>(m_N);
    m_gpep = zeros<vec>(m_N);

    m_sd1 = zeros<vec>(m_N);
    m_sd2 = zeros<vec>(m_N);
    m_stn = zeros<vec>(m_N);
    m_gpi = zeros<vec>(m_N);
    m_gpe = zeros<vec>(m_N);
}


////////////////////////////////////////////////////////////////////////////

void BasalGanglia::save()
{
    // save storage data to file
    
    ofstream fout( (m_filename+ "_save").c_str() );

    for(int r=0; r<m_N*m_NUMCOMP; r++ )
    {
        for(int c=0; c<m_STIME; c++)
        {
            fout << m_data(r,c) << " ";
        }
        fout << endl;
    }
}

////////////////////////////////////////////////////////////////////////////


void BasalGanglia::step(vec inp, vec crx, vec da)
{

    // A step of numerical integration. Each layer is updated using euler integration:
    //
    // the leaky integrator        tau*dy/dt = -y + input + W*f(y) 
    // becomes                     y(t) = y(t-1) + (dt/tau) * (-y(t-1) + input +  W*f(y))


    double h = m_DT/m_TAU;

    // striatum D1
    m_sd1p += h*( - m_sd1p
            + (m_SD1_BL+m_SD1_DA*da) % (    // dopamine multiplies: (a + Da) * input 
                + m_INP2SD1_W*inp 
                + m_CRX2SD1_W*crx ) );

    m_sd1 = outfun(m_sd1p, m_TH, m_AMP);

    // striatum D2
    m_sd2p += h*( - m_sd2p
            + ( 1.0/(m_SD2_BL+m_SD2_DA*da) ) % (    // dopamine divides: (1/(a + Da)) * input
                + m_INP2SD2_W*inp
                + m_CRX2SD2_W*crx ) );
    m_sd2 = outfun(m_sd2p, m_TH, m_AMP);

    // STN
    m_stnp += h*( - m_stnp
            + m_INP2STN_W*inp
            + m_CRX2STN_W*crx
            - m_GPE2STN_W*m_gpe );
    m_stn = outfun(m_stnp, m_TH, m_AMP);

    // GPe
    m_gpep += h*( - m_gpep
            + m_GPE_BL
            - m_SD22GPE_W*m_sd2
            + (ones<mat>(m_N,m_N)*m_STN2GPE_W) * m_stn );

    m_gpe = outfun(m_gpep, m_TH, m_AMP);
    
    // GPi
    m_gpip += h*( - m_gpip
            + m_GPI_BL
            - m_SD12GPI_W*m_sd1
            + (ones<mat>(m_N,m_N)*m_STN2GPI_W) * m_stn 
            - m_GPE2GPI_W*m_gpe);

    m_gpi = outfun(m_gpip, m_TH, m_AMP);

}

////////////////////////////////////////////////////////////////////////////


void BasalGanglia::store(int t)
{
    int c = 0;
    
    //store current activations
    m_data(span(m_N*c,m_N*(c+1)-1), t) = m_sd1; c++; 
    m_data(span(m_N*c,m_N*(c+1)-1), t) = m_sd2; c++;
    m_data(span(m_N*c,m_N*(c+1)-1), t) = m_stn; c++;
    m_data(span(m_N*c,m_N*(c+1)-1), t) = m_gpi; c++;
    m_data(span(m_N*c,m_N*(c+1)-1), t) = m_gpe;
}

////////////////////////////////////////////////////////////////////////////

