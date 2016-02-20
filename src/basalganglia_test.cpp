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

#include <basalganglia.h>
#include <layer.h>


int main()
{

    srand(time(0));

    try
    {
        // initialize objects
        BasalGanglia bg("basalganglia");
        Layer tha("thalamus");
        Layer crx("cortex");
        
       
        // define constants

        int STIME = bg.m_STIME;
        int N = bg.m_N;

        double INP_GAP = .8;
        double INP_STDERR = .1;
        double INP_BLERR = .2;
        double CRX_STDERR = .4;
        double CRX_BLERR = .01;

        double DA_MIN = .5;
        double DA_MAX = 1;

        // define weights
        double GPI2THA_W = 3;
        double THA2CRX_W = 4;
        double CRX2THA_W = .8;

        // create the timeseries of inputs (only noise)
        mat inp = INP_STDERR*randn<mat>(N,STIME)+ INP_BLERR;
        
        // create the timeseries of dopamine activity (all values to da_min)
        mat da = DA_MIN*ones<vec>(N)*ones<vec>(STIME).t();

        // add scheduling variations to input and dopamine:
        // the first half of the simulation is  divided in three intervals,
        //
        for(int t=0; t<STIME; t++)
        {         
           
            // add a further measure to the first channel in the 
            // first interval
            if ( (STIME*(1/12.-1/24.))  < (t%(STIME/2))  and 
                    (t%(STIME/2)) < (STIME*(1/12.+1/24.)) )
                inp(0,t) += INP_GAP;

            // add a further measure to the second channel in the 
            // second interval 
            if ( (STIME*(3/12.-1/24.)) < (t%(STIME/2)) and 
                    (t%(STIME/2)) < (STIME*(3/12.+1/24.)) )  
                inp(1,t) += INP_GAP;
          
            // add a further measure to the third channel in the 
            // third interval 
            if ( (STIME*(5/12.-1/24.)) < (t%(STIME/2)) and
                    (t%(STIME/2)) < (STIME*(5/12.+1/24.)) )
                inp(2,t) += INP_GAP;


            // add a further measure to dopamine for each interval
            if ( (STIME/12) < (t%(STIME/6)) and 
                    (t%(STIME/6)) < (STIME*(5/12.)) )  
                da(span::all,t).fill(DA_MAX);

            // dopamine stays high during the second half of
            // the simulation
            if ( t > (STIME/2 + STIME/12) ) 
                da(span::all,t).fill(DA_MAX);

        }

        // spreading
        for(int t=0; t<STIME; t++)
        {
            
            bg.step(
                    inp.col(t),
                    crx.m_output,
                    da.col(t)
                    );
            
            tha.step( 
                    CRX2THA_W*crx.m_output - 
                    GPI2THA_W*bg.m_gpi
                    );
           
            crx.step( 
                    THA2CRX_W*tha.m_output + 
                    CRX_BLERR + 
                    CRX_STDERR*randn<vec>(N) 
                    );

            bg.store(t);
            tha.store(t);
            crx.store(t);

        }

        // print data of the  simulation as a 
        // (STRD1+STRD2+STN+GPi+GPe+THA+CRX)xSTIME matrix
        join_vert(
                join_vert(
                    bg.m_data,
                    tha.m_data),
                crx.m_data).raw_print();

    }
    catch(parameter_exception &e)
    {
        e.what();
        return 1;
    }

    return 0;
}
