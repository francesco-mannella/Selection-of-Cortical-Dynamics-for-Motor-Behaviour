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

#ifndef PARAMETRIZABLE_H
#define PARAMETRIZABLE_H

#include <memory>
#include <stdexcept>

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>
#include <armadillo>

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>


////////////////////////////////////////////////////////////////////////////////
/// Exceptions /////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Base class of Parametrizable exceptions
 */
class parameter_exception : public std::runtime_error 
{
    public:    
        parameter_exception(const std::string &msg) 
            : std::runtime_error(msg) { }
};

/**
 * @brief parameter-not-found exception
 *
 * Thrown when a parameter name is not found 
 * in the parameters file 
 */
class parameter_name_exception : public parameter_exception 
{
    public:    
        parameter_name_exception(const std::string &name)
            : parameter_exception(std::string("Parameter '") 
                    +name + "' not found!") {}
};

/**
 * @brief file-not-found exception
 *
 * Thrown when the parameters file is not found 
 */
class parameter_file_exception : public parameter_exception 
{
    public:    
        parameter_file_exception(const std::string &file)
            : parameter_exception(std::string("File '") 
                    +file + "' not found!") {}
};

/**
 * @brief dir-not-found exception
 *
 * Thrown when the m_parameters directory cannot be created  
 */
class parameter_dir_exception : public parameter_exception 
{
    public:    
        parameter_dir_exception(const std::string &dir)
            : parameter_exception(std::string("Could not create directory '") 
                    +dir +"'!") {}
};


////////////////////////////////////////////////////////////////////////////////
// iO operators
// VECTOR 

// overloading of operator "<<" for std::vector 
template<typename T>
std::ostream& operator<< (
        std::ostream& stream,  
        std::vector<T> &v)
{
    static_assert(std::is_arithmetic<T>::value,"T must be a number");
    
    stream << "[";
    for(int x=0; x<v.size()-1;x++ )
        stream << v[x] << ",";
    stream << v[v.size()-1] << "]";

    return stream;
}

////////////////////////////////////////////////////////////////////////////////

template<typename T> std::istream& operator>> (
        std::istream& stream,  
        std::vector<T> &v)
{
    static_assert(std::is_arithmetic<T>::value,"T must be a number");


    std::string r;
    std::getline(stream,r);

    size_t end = r.find("]");
    size_t begin = r.rfind("[",end);

    if(begin != std::string::npos and end != std::string::npos )
    { 
        T tmp;
        v.clear();
        r = std::string(r.begin()+begin+1,r.begin()+end); 
        std::string e;
        std::stringstream substream(r);
        while(std::getline(substream,e,',')) 
        {
            std::stringstream(e) >> tmp;
            v.push_back(tmp);
        }
    }

    return stream;
}

////////////////////////////////////////////////////////////////////////////////

inline std::istream& operator>> (
        std::istream& stream,  
        arma::vec &v)
{

    std::string r;
    std::getline(stream,r);

    size_t end = r.find("]");
    size_t begin = r.rfind("[",end);

    if(begin != std::string::npos and end != std::string::npos )
    { 
       
        double tmp;
        std::vector<double> vt;
        r = std::string(r.begin()+begin+1,r.begin()+end); 
        std::string e;
        std::stringstream substream(r);
        while(std::getline(substream,e,',')) 
        {
            std::stringstream(e) >> tmp;
            vt.push_back(tmp);
        }
        v = vt;

    }

    return stream;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Base of Parameter class for factorization
 *
 * Parameter objects of different types inherit from 
 * a common base so that they can be collected in 
 * a single container.
 */
class BaseParameter 
{
    public:

        virtual void update(const std::vector<std::string> &lines ) = 0;
        virtual void save(std::ostream &out) = 0;
    

};

////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Generic Parameter class 
 *
 * A generic definition so that Parameter objects can be linked to variables 
 * of different types. The only constraint is that a type is recognized by 
 * the std::ostream operator '<<' and the std::istream operator '>>'.
 */
template <typename T> class Parameter : public BaseParameter 
{

    public:

        /**
         * @param   nname       Name of the parameter 
         * @param   value       Reference to the variable to be synchronized 
         * @param   separator   value - label separator in the parameter file 
         */
        Parameter<T>(const std::string &name,  T &value, 
                const std::string &separator ="-------") 
            : m_name(name),m_value(value),m_separator(separator) {}


        // Avoid copying or assignement ////////////////////////////////////////
        
        Parameter<T>() = delete;
        Parameter<T>(const Parameter<T>& ) = delete;
        Parameter<T>& operator=(const Parameter<T>& ) = delete; 

        ////////////////////////////////////////////////////////////////////////

        /**
         * Update the values of the m_parameters 
         *
         * @param   lines   the lines of the parameter file
         *                  parsed to search values
         */
        virtual void update( const std::vector<std::string> &lines );
        
        /**
         * Save the values of the m_parameters 
         *
         * @param   out     Storage stream
         */ 
        virtual void save(std::ostream &out);

    private :

        const std::string m_name;    /*!< Name of the parameter */
        const std::string m_separator;  /*!< separator */
        T &m_value; /*!< reference to actual value */

};



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


/**
 * @brief Allows to synchronize the values of variables with 
 * those of a parameter file
 *  
 * A Parametrizable object is able to manage the synchronization 
 * between variables in the program and a parameter file. 
 * 
 * usage:
 *
 * \code
 *  
 *  using namespace std;
 *
 *  ...
 *
 *  // Some variables
 *  int first_variable;
 *  double second_variable;
 *  std::string third_variable;
 *
 *  ...
 *
 *  // Create a Parametrizable object
 *  Parametrizable params("a_parameter_file");
 *
 *  ...
 *
 *  // Link variables to the object
 *  params.addParameter("first_variable_LABEL", first_variable);
 *  params.addParameter("second_variable_LABEL", second_variable);
 *  params.addParameter("third_variable_LABEL", third_variable);
 *
 *  ...
 *
 *  try
 *  {
 *      params.loadParameters(); // in case of success the values 
 *                               // first_variable, second_variable,
 *                               // and third_variable are updated
 *                               // according to the values stored 
 *                               // inthe file params.m_filedir
 *  }
 *  catch(parameter_name_exception &e)    // check for incongruence:
 *  {                                     // parameter file exists but 
 *      cout << e.what() << endl;         // not all linked variables 
 *      exit(1);                          // have a stored value
 *  }
 *  catch(parameter_file_exception &e)   // check for existence:
 *  {                                    // the file do not exists so
 *                                       // we can create one
 *      // Initiale variables
 *      first_variable = 1;             
 *      second_variable = 2.0;
 *      third_variable = "three";
 *      
 *      // Create m_parameters file
 *      // with the current values
 *      params.saveParameters();
 *  }
 *
 *
 * \endcode
 */
class Parametrizable
{

    public :

        const std::string m_filedir;    /*!< Directory of the 
                                            parameter file */
        const std::string m_filename;   /*!< Name of the 
                                            parameter file */

        /**
         * @param   filename  Parameter file name 
         * @param   filedir  Parameter file directory 
         */
        Parametrizable(const std::string &filename, 
                const std::string &filedir = "parameters" ) 
            : m_filename(
                    boost::filesystem::path(std::string(filedir)+"/"+filename).string() ),
            m_filedir(filedir)
        { 
            // Throw an exception if could not create parameter file directory
            if (not boost::filesystem::exists(m_filedir))
                if(not boost::filesystem::create_directory(m_filedir))
                    throw parameter_dir_exception(m_filedir);
        }
       
        // Avoid copying or assignement ///////////////////////////////////////
        
        Parametrizable() = delete;
        Parametrizable(const Parametrizable& ) = delete;
        Parametrizable& operator=(const Parametrizable& ) = delete;   

        ////////////////////////////////////////////////////////////////////////

        /**
         * Add a Parameter object that is linked to an exernal variable
         *
         * @param   name    Name of the parameter
         * @param   value   Reference to the variable to be synchronized 
         */
        template <typename T> void addParameter(
                const std::string &name,
                T &value);  

        /**
         * Load values from the parameter file and update 
         * the linked variables
         */
        void loadParameters();

        /**
         * Save values to the parameter file and
         */
        void saveParameters();
    

    private:


        /** Container of created Parameter objects */
        std::vector<std::unique_ptr<BaseParameter>> m_parameters;

};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Definitions of Parameter<T> methods

template <typename T> void inline Parameter<T>::update( 
        const std::vector<std::string> &lines ) 
{

    // pattern for the label  
    std::string pattern = 
        std::string("\\s*") + m_separator + "\\s*\\<" + 
        m_name + "\\>\\s*$";

    boost::regex re(  pattern );   

    // extract value from lines
    bool parameter_found = false;
    for( auto &line : lines  ) 
    {
        if ( boost::regex_search( line, re ) ) 
        {  
            std::string rep =  
                boost::regex_replace( line , re, "" );
            std::stringstream buffer(rep);
            buffer >> m_value;

            parameter_found = true;
            break;
        }
    }

    if (not parameter_found)
        throw parameter_name_exception(m_name);

}

////////////////////////////////////////////////////////////////////////////////

template <typename T> void inline Parameter<T>::save(std::ostream &out) 
{
    out << m_value << " " << m_separator << " " << m_name << std::endl;
};



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Definitions of Parametrizable methods

template <typename T> void inline Parametrizable::addParameter(
        const std::string &name, T &value)
{
    m_parameters.push_back( std::unique_ptr<BaseParameter>(
                new Parameter<T>(name,value)) );
}

////////////////////////////////////////////////////////////////////////////////

void inline Parametrizable::loadParameters( )
{
    
    // Open parameter file
    std::ifstream par_file(m_filename.c_str()); 
    if(not par_file.is_open())
        throw parameter_file_exception(m_filename);

    // Load lines from the file
    std::vector<std::string> lines;
    std::string line("");
    while( std::getline( par_file, line ) )
        lines.push_back( line );
    lines.push_back( line );

    // Update values
    for(auto &parameter: m_parameters )
        parameter->update(lines);
}

////////////////////////////////////////////////////////////////////////////////

void inline Parametrizable::saveParameters()
{
    // Open parameter file
    std::ofstream par_file(m_filename.c_str()); 
    
    // Save values
    for(auto &parameter: m_parameters )
        parameter->save(par_file);
}


#endif //PARAMETRIZABLE_H


