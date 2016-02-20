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


int main()
{

    srand(time(0));

    Reservoir r("reservoir");
    r.init();

    int STIME = r.m_STIME;
    int N = r.m_N;

    vec rows = zeros(N); 
    mat inp = zeros(N,STIME);
    inp.col(0) = randn(N);

    ofstream data_o("activation_data");
    ofstream w_o("weights");

    for(int t=0; t<STIME; t++)
    {
    
        r.step(inp.col(t));
        r.m_output.t().raw_print(data_o);
    }
   
    r.m_w.raw_print(w_o);

    return 0;
}
