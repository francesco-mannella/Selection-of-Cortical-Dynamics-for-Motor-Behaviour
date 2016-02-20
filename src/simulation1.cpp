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

#include <simulator.h>

// Implements simulation in Fig. 7

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
// MAIN //////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  
    int SEED = 1;   // Seed - 1 if not updated by command-line arg
    bool ONLINE = false;   // Flag for online learning 

    //////////////////////////////////////////////////////////////////////////////////////////
    // MANAGE MAIN ARGUMENTS ///////////////////////////////////////////////////////////////// 
    //////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    enum OPTIONS { O_HELP=500, O_ONLINE, O_SEED };    // Indices of command line options
    
    map<string,int> option;    // Maps command-line strings to indices
    option["-h"] = O_HELP;
    option["--help"] = O_HELP;
    option["-o"] = O_ONLINE;
    option["--online"] = O_ONLINE;
    option["-s"] = O_SEED;
    option["--seed"] = O_SEED;


    // Help string to show in case of exit
    string help =
        " usage csntc_sim [options] \n\n"
        "   -h,   --help            show this help \n"
        "   -o,   --online          set online learning \n"
        "   -s N, --seed N          seed \n"
        "\n";

    // Iterate over args
    for(int x=1; x<argc; x++)
    {
        // Find the corresponding option
        if ( option.find(argv[x]) == option.end() )
        {      
            // Control for consistency (for options that require extra arguments) 
            if( not ( option[argv[x-1]] == O_SEED and is_number(argv[x]) ) )       
            {
                cout << endl << "wrong parameters" << endl << endl;
                cout << help; exit(0);
            } 
        }

        // Do the proper action for each option
        switch(option[argv[x]]) 
        {
            case O_HELP :     // Print help string and exit
                cout << help; 
                exit(0); 
                break;
   
            case O_ONLINE :     // Switch online flag to true
                ONLINE = true;
                break;

            case O_SEED :     // Update seed
                try 
                {
                    str2type( argv[x+1], SEED); 
                } 
                catch (exception &ex) 
                { 
                    cout << "wrong argument" << endl << endl;
                    cout << help; 
                    exit(1);
                }
                break;
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    // Print updated seed

    cout << endl;
    cout << "ONLINE: " << (ONLINE?"TRUE":"FALSE") << endl;
    cout << "SEED: " << SEED << endl;
    cout << endl;

    //////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////
    

    // Run offline learning
    try 
    {
        srand(SEED);
        Simulator sim("simulation1",SEED);
        srand(SEED);    // Reset seed to avoid differences due to the constructor
        if (ONLINE)
            sim.run_online();
        else
            sim.run_offline();
    }
    catch(parameter_exception &e)
    {
        e.what();
        exit(1);
    }


}

