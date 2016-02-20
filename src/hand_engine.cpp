#include "hand_engine.h"

using namespace cens;


// Then define the initObjects function
// in which we initialize the rigid body
// as a box
void Engine::initObjects() {
    
    // Simply import all objects form a *.bullet file
    
    m_eRobots["hand"]=loadBulletFile("hand","hand.bullet");

    // Update collision table

    CENSNonCollidingTable table;
    
    table[CENSNCPair(  m_eBodies["Palm"], m_eBodies["Wrist"] ) ] = false;
    table[CENSNCPair(  m_eBodies["Palm"], m_eBodies["PT"] ) ] = false;
    table[CENSNCPair(  m_eBodies["PT"], m_eBodies["Thumb1"] ) ] = false;
    table[CENSNCPair(  m_eBodies["Palm"], m_eBodies["Thumb1"] ) ] = false;

    table[CENSNCPair(  m_eBodies["Palm"], m_eBodies["PI"] ) ] = false;
    table[CENSNCPair(  m_eBodies["PI"], m_eBodies["Index1"] ) ] = false;
    table[CENSNCPair(  m_eBodies["Palm"], m_eBodies["Index1"] ) ] = false;

    table[CENSNCPair(  m_eBodies["Palm"], m_eBodies["PM"] ) ] = false;
    table[CENSNCPair(  m_eBodies["PM"], m_eBodies["Miiddle1"] ) ] = false;
    table[CENSNCPair(  m_eBodies["Palm"], m_eBodies["Miiddle1"] ) ] = false;

    table[CENSNCPair(  m_eBodies["Palm"], m_eBodies["PR"] ) ] = false;
    table[CENSNCPair(  m_eBodies["PR"], m_eBodies["Ring1"] ) ] = false;
    table[CENSNCPair(  m_eBodies["Palm"], m_eBodies["Ring1"] ) ] = false;

    table[CENSNCPair(  m_eBodies["Palm"], m_eBodies["PP"] ) ] = false;
    table[CENSNCPair(  m_eBodies["PP"], m_eBodies["Pinky1"] ) ] = false;
    table[CENSNCPair(  m_eBodies["Palm"], m_eBodies["Pinky1"] ) ] = false;

    setCollisionFilter(table);
    toggleTexture();
  
    // set initial color of all objects
    
    CENSSerializedRobot::Bodies bds = m_eRobots["hand"]->bodies;
    CENSSerializedRobot::Bodies::iterator it = bds.begin();
    for( ; it != bds.end(); ++it) 
        m_eShapes[ (*it).second ] ->setColor(Vector3f(.2,.2,.2));
    
    // initialy disable all bodies (to avoid body-sleep)
    m_eRobots["hand"]->disableBodies();       
    m_eRobots["hand"]->enableBodies();       

    // allocate the controller
    int hinge_number = m_eRobots["hand"]->hinges.size()-1; 
    controller = shared_ptr<Controller>( new Controller(seed, hinge_number) );

    savestream.open("save");
}


// Finally, define the step function
// in which we apply a force to the box
// for the first 10 steps of the simulation
void Engine::step( int timestep ) 
{
    // setting environment variables 
    
    // counters 
    static int  t = 0;     // within-trial timestep count
    static int  tt = 0;    // trial count   

    //static const int after_trial_interval = 1;
    static const int trial_interval = 1000;

    // setting controller variables 
    const double impulse = .8;
    const int learn_trials = 3*10;

    ///////////////////////////////////////////////////////////////////////////// 
    ///////////////////////////////////////////////////////////////////////////// 
    ///////////////////////////////////////////////////////////////////////////// 
    
    // read current hinge angles

    vec real_patterns;

    real_patterns <<          
        m_eRobots["hand"]->hinges["PT_to_Thumb1"].hinge->getHingeAngle()<<
        m_eRobots["hand"]->hinges["Thumb1_to_Thumb2"].hinge->getHingeAngle()<<
        m_eRobots["hand"]->hinges["Thumb2_to_Thumb3"].hinge->getHingeAngle()<<
        m_eRobots["hand"]->hinges["PI_to_Index1"].hinge->getHingeAngle()<<
        m_eRobots["hand"]->hinges["Index1_to_Index2"].hinge->getHingeAngle()<<
        m_eRobots["hand"]->hinges["Index2_to_Index3"].hinge->getHingeAngle()<<
        m_eRobots["hand"]->hinges["PM_to_Middle1"].hinge->getHingeAngle()<<
        m_eRobots["hand"]->hinges["Middle1_to_Middle2"].hinge->getHingeAngle()<<
        m_eRobots["hand"]->hinges["Middle2_to_Middle3"].hinge->getHingeAngle()<<
        m_eRobots["hand"]->hinges["PR_to_Ring1"].hinge->getHingeAngle()<<
        m_eRobots["hand"]->hinges["Ring1_to_Ring2"].hinge->getHingeAngle()<<
        m_eRobots["hand"]->hinges["Ring2_to_Ring3"].hinge->getHingeAngle()<<
        m_eRobots["hand"]->hinges["PP_to_Pinky1"].hinge->getHingeAngle()<<
        m_eRobots["hand"]->hinges["Pinky1_to_Pinky2"].hinge->getHingeAngle()<<
        m_eRobots["hand"]->hinges["Pinky2_to_Pinky3"].hinge->getHingeAngle();

    real_patterns.t().raw_print(savestream);

    ///////////////////////////////////////////////////////////////////////////// 
    ///////////////////////////////////////////////////////////////////////////// 
    ///////////////////////////////////////////////////////////////////////////// 
    
    // set dopamine and switch learning
    double dopamine  =
        (t>(trial_interval*(1/8.)) and t< trial_interval*(7./8.))? 1.0 : 0.0; 
    bool learn = false; if (tt < learn_trials) learn = true;   
    bool learn_int = false;
    if (learn == true and t >(trial_interval*(1/2.)) and t<(trial_interval*(3/4.)))
        learn_int = true;
    learn_int = learn_int and learn;

    // change colors
    CENSSerializedRobot::Bodies bds = m_eRobots["hand"]->bodies;
    CENSSerializedRobot::Bodies::iterator it = bds.begin();
    if (dopamine == 0.0)
        for( ; it != bds.end(); ++it) 
            m_eShapes[ (*it).second ] ->setColor(Vector3f(.2,.2,.2));
    else if (dopamine == 1.0)
        if (learn_int == false)
            for( ; it !=bds.end(); ++it) 
                m_eShapes[ (*it).second ] ->setColor(Vector3f(.5,.5,.5));
        else
            for( ; it !=bds.end(); ++it) 
                m_eShapes[ (*it).second ] ->setColor(Vector3f(1,.2,.2));


    vec readout = (M_PI/2.)*controller->step(t, tt,trial_interval, dopamine, learn);

    ///////////////////////////////////////////////////////////////////////////// 
    ///////////////////////////////////////////////////////////////////////////// 
    ///////////////////////////////////////////////////////////////////////////// 

    // real motor commands
    int i=0;
    m_eRobots["hand"]->move("PT_to_Thumb1", readout(i++), impulse);
    m_eRobots["hand"]->move("Thumb1_to_Thumb2", readout(i++), impulse);
    m_eRobots["hand"]->move("Thumb2_to_Thumb3", readout(i++), impulse);
    m_eRobots["hand"]->move("PI_to_Index1", readout(i++), impulse);
    m_eRobots["hand"]->move("Index1_to_Index2", readout(i++), impulse);
    m_eRobots["hand"]->move("Index2_to_Index3", readout(i++), impulse);
    m_eRobots["hand"]->move("PM_to_Middle1", readout(i++), impulse);
    m_eRobots["hand"]->move("Middle1_to_Middle2", readout(i++), impulse);
    m_eRobots["hand"]->move("Middle2_to_Middle3", readout(i++), impulse);
    m_eRobots["hand"]->move("PR_to_Ring1", readout(i++), impulse);
    m_eRobots["hand"]->move("Ring1_to_Ring2", readout(i++), impulse);
    m_eRobots["hand"]->move("Ring2_to_Ring3", readout(i++), impulse);
    m_eRobots["hand"]->move("PP_to_Pinky1", readout(i++), impulse);
    m_eRobots["hand"]->move("Pinky1_to_Pinky2", readout(i++), impulse);
    m_eRobots["hand"]->move("Pinky2_to_Pinky3", readout(i++), impulse);
    m_eRobots["hand"]->move("Palm_to_PT",readout(i++), impulse);
    m_eRobots["hand"]->move("Palm_to_PI",readout(i++), impulse);
    m_eRobots["hand"]->move("Palm_to_PM",readout(i++), impulse);
    m_eRobots["hand"]->move("Palm_to_PR",readout(i++), impulse);
    m_eRobots["hand"]->move("Palm_to_Pp",readout(i++), impulse);

    // uncomment to save frames as files
    //if ((trial_interval*tt + t) >0) saveCameraPixelMap(0, trial_interval*tt + t,"png");
    
    t++;
    if(t>trial_interval) 
    {
        tt++;
        t = 0;
    }


    // This call should always be here. 
    // Without it bullet and graphics steps 
    // are not called, and the loop of 
    // simulation doesn't go on.
    CENSSerializedEngine::step(timestep);

}

