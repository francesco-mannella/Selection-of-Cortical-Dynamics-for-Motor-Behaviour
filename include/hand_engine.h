#ifndef ENGINE_H
#define ENGINE_H

#include <map>
#include <fstream>
using namespace std;


#include "CENS/cens_serialized_engine.h"
#include "hand_controller.h"
using namespace cens;



// First, declare a CENSEngine derived class
class Engine : public CENSSerializedEngine {   

    public:   
        
        Engine( unsigned int _seed):
            CENSSerializedEngine(),seed(_seed){}
        
        void initObjects();
        void step( int timestep);

        unsigned int seed;
       
        // smart pointer to controller
        shared_ptr<Controller> controller;

        // data stream
        ofstream savestream;
        
};

#endif //ENGINE_H

