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

#include "layer.h"

//////////////////////////////////////////////////
//////////////////////////////////////////////////


Layer::Layer(std::string name) :       
    Parametrizable(name)
{
        addParameter("STIME",m_STIME);
        addParameter("DT",m_DT);
        addParameter("TAU",m_TAU);
        addParameter("TH",m_TH);
        addParameter("AMP",m_AMP);
        addParameter("N",m_N);
        addParameter("BL",m_BL);

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
        m_BL = 0.0;

        saveParameters();
        
    }


    reset();
    
    m_data = zeros<mat>(m_N,m_STIME);
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////

void Layer::step(vec inp)
{
    double h = m_DT/m_TAU;

    m_u += h*( - m_u + m_BL + inp );
    m_v = outfun(m_u, m_TH, m_AMP);

    m_output = m_v;

}

//////////////////////////////////////////////////
//////////////////////////////////////////////////

void Layer::store(int t)
{
    m_data(span::all, t) = m_output;

}

//////////////////////////////////////////////////
//////////////////////////////////////////////////

void Layer::reset()
{
    m_u = zeros<vec>(m_N);
    m_v = zeros<vec>(m_N);
    m_output = zeros<vec>(m_N);
}
