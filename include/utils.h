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

#ifndef UTILS_H
#define UTILS_H

#include <type_traits>
#include <iostream>
#include <armadillo>
#include <exception>
#include <sstream>
#include <map>


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UTILITIES ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using namespace std;
using namespace arma;


////////////////////////////////////////////////////////////////////////////////

/** Finds if a type is an armadillo array - metafunction  */

template< typename, typename = void >
struct is_armadillo
  : std::false_type {};


template< typename T >
struct is_armadillo< T, 
    typename std::enable_if<std::is_member_function_pointer<decltype(&T::t)>::value>::type >
    : std::true_type {};




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


vec outfun(vec x, double th = 0, double amp = 1);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


vec sinusoid_sequence( 
        double steps=100,  
        double period = 1,
        double phase = 0,
        double center = 0,
        double amplitude = 1 );

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * NRMSE - Normalized Root Mean Square Error
 *
 * @param   f       Predicted values
 * @param   y       Expected values
 */
double mse(mat f, mat y);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Identify a string as an integer or floating-point number
 *
 * @param   s       A string 
 */
bool is_number(const std::string& s);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Tests if the competition in the basal ganglia is optimal for all the three trials in a session 
 *
 * @param   tha     Matrix of activations in the thalamus
 * @param   rng1    Range of activity of the first trial
 * @param   rng2    Range of activity of the second trial
 * @param   rng3    Range of activity of the third trial
 */
bool test_bg_competition(mat tha, umat rngs );

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * convert index to a unitary index vector 
 */
uvec uidx(int index);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * Converts from string to value 
 *
 * @param   str     A string containing a value of type T
 * @param   el      A reference to a variable where to put the converted element
 */
template<typename T>
void str2type(const string &str, T &el)
{
    stringstream buf;
    
    // Allows buf to throw exceptions
    buf.exceptions(ios::failbit | ios::badbit);
    
    buf.str(str);
    buf >> el;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Extend division reminder to vectors
 *
 * @param   a       Dividend 
 * @param   n       Divisor
 */
template<typename T>
T mod(T a, int n)
{
    return a - floor(a/n)*n;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef NOLOG

// empty one-argument version
template<typename T>
inline void print_log(T&& t) {  }

// empty one-argument-no-endline version
template<typename T>
inline void print_log_not_ended(T&& t ) { }

#else

// one-argument version
template<typename T>
inline void print_log(T&& t ) { cout << t << endl;  }

// one-argument-no-endline version
template<typename T>
inline void print_log_not_ended(T&& t ) { cout << t; }

#endif

// termination version
inline void print_log() {  cout << endl;  }

// variadic version
    template<typename Arg, typename... Args>
inline void print_log(Arg&& arg1, Args&&... args)
{
    print_log_not_ended(arg1);
    print_log(args...);
    if( sizeof...(args) < 1  )
        print_log();
}

// reminder (custom version)
inline double remind(double a,double b)
{
    return (a/b) -((int)a/(int)b);
}

////////////////////////////////////////////////////////////////////////////////

/**  Element-wise maximum  */

template< typename T1, typename T2,
    typename std::enable_if<std::is_arithmetic<T2>::value>::type* = nullptr>
T2 maximum(T1 th, T2 v)
{
    static_assert( std::is_arithmetic<T1>::value == true, 
            "First argument is not arithmentic" );
    return v*(v>=th) + th*(v<th);
}

template< typename T1, typename T2,
    typename std::enable_if<is_armadillo<T2>::value>::type* = nullptr>
T2 maximum(T1 th, T2 v)
{
    static_assert( std::is_arithmetic<T1>::value == true, 
            "First argument is not arithmentic" );
    return v%(v>=th) + th*(v<th);
}


////////////////////////////////////////////////////////////////////////////////

/**  Element-wise minimum  */

template< typename T1, typename T2,
    typename std::enable_if<std::is_arithmetic<T2>::value>::type* = nullptr>
T2 minimum(T1 th, T2 v)
{
    static_assert( std::is_arithmetic<T1>::value == true, 
            "First argument is not arithmentic" );
    return v*(v<=th) + th*(v>th);
}

template< typename T1, typename T2,
    typename std::enable_if<is_armadillo<T2>::value>::type* = nullptr>
T2 minimum(T1 th, T2 v)
{
    static_assert( std::is_arithmetic<T1>::value == true, 
            "First argument is not arithmentic" );
    return v%(v<=th) + th*(v>th);
}

////////////////////////////////////////////////////////////////////////////////

/**  Element-wise bound  */

    template< typename T1, typename T2>
T2 bound(T1 min, T1 max, T2 v)
{
    static_assert( std::is_arithmetic<T1>::value == true, 
            "First argument is not arithmentic" );
    return maximum(min,minimum(max,v));
}



#endif //UTILS_H

