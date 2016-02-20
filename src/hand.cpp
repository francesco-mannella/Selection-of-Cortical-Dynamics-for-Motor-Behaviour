#include <fenv.h>
#include "hand_engine.h"

using namespace cens;

int main(int argc, char** argv) 
{
    
    // set seed to current time
    unsigned int rnd = time(0);

    // change seed to the one given by user
    if(argc == 2 ) stringstream(argv[1]) >> rnd;
    
    // Trap floating-point exceptions
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
   
    // initialize random number generation
    srand(rnd);

    // Create and initialize engine
    Engine engine(rnd);
    engine.init(argc,argv);

    // start simulation loop
    engine.run();

    return 0;
}

