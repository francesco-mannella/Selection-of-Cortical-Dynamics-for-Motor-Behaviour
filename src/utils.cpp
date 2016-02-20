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

#include <utils.h>

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vec outfun(vec x, double th, double amp) 
{
    vec y = tanh(amp*(x-th));
    return y;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vec sinusoid_sequence( 
        double steps,  
        double period,
        double phase,
        double center,
        double amplitude ) 
{
    vec x = linspace(0,steps-1,steps);
    return amplitude*cos( period*2*datum::pi*(x/steps) + phase ) + center;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double mse(mat f, mat y) 
{
    return sqrt( accu( pow((f-y),2) )  /f.n_elem ) / (y.max()-y.min());
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();   
    while (it != s.end() && (std::isdigit(*it) or *it=='.' ) ) ++it;     // Iterate through the string characters, stops if not [.0123456789] 
    return not s.empty() && (it == s.end());    // Return true if it is the terminator of a non-epmty string 
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool test_bg_competition(mat tha, umat rngs )
{
    bool test_channel = true;
    int n_rngs = rngs.n_cols;
    uvec test_rngs = zeros<uvec>(n_rngs);
    double test_threshold = 0.1;

    for(int r=0; r<tha.n_rows; r++)    // For each channel
    { 
     
        vec mm = zeros<vec>(n_rngs);
        for(int c=0; c<n_rngs; c++)
        {
            // Compute the means of activities during the three different ranges
            double m = mean(rowvec(tha.row(r))(rngs.col(c)));

            // Test if no more than a channel wins during each trial
            test_rngs(c) = (m>test_threshold and (not test_rngs(c))) or (m<test_threshold and test_rngs(c) );

            mm[c] = m;
        }

        // Test if each channel wins no more than a trial over three in a session
        test_channel = test_channel  and ( accu(mm>test_threshold) == 1 );

    }

    // All conditions must be true (1/3 channels per each trial and 1/3 wins per channel in the session)
    return all(test_rngs > 0 ) and test_channel;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

uvec uidx(int index)
{
    return ones<uvec>(1)*index;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

