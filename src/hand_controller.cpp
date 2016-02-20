#include "hand_controller.h"

Controller::Controller(int s, int n)
    :
        n_readout(n),
        csntc(new CSNTC("csntc",s))
{
    readout = zeros<vec>(n_readout);
    readout_pot = zeros<vec>(n_readout);
    
    inp = zeros<vec>(csntc->m_N_CRX);
    inp_prev = zeros<vec>(csntc->m_N_CRX);
    
    err = zeros<vec>(n_readout);;
    err_prev = zeros<vec>(n_readout);
    
    activity = zeros<vec>(csntc->m_N_CRX);;
    activity_prev = zeros<vec>(csntc->m_N_CRX);

    W_readout = randn<mat>( n_readout, csntc->m_N_CRX );

    W_env2CRX = zeros(csntc->m_N_CRX);
    W_env2CRX(  span(csntc->m_N_CRX - csntc->m_THA2CRX_SPREAD, csntc->m_N_CRX-1 ) ) = 1;

    trialstream.open("trial");
}


const vec &Controller::step(int t, int tc, int ti, double d, bool learn)
{
   
    static double dt = 0.001;
    static double decay = 0.005;
    static double delta = dt/decay;
    static mat target = mat(
            ".95 .95 .95 .00 .00 .00 .95 .95 .95 .95 .95 .95 .95 .95 .95 .00 .00 .00 .00 .00;"
            ".95 .95 .95 .00 .00 .00 .00 .00 .00 .95 .95 .95 .95 .95 .95 .00 .00 .00 .00 .00;"
            ".95 .95 .95 .00 .00 .00 .00 .00 .00 .00 .00 .00 .95 .95 .95 .00 .00 .00 .00 .00" );
    
    static uvec testtimeline = uvec("0 1 2");
 
    
    if(tc%trials==0 and t==0) testtimeline = shuffle(testtimeline);
    
    //
    // SPREAD
    //
    trialstream << testtimeline[tc%trials] << endl;
    vec bg = zeros(trials);  bg(testtimeline[tc%trials]) = 1;
    vec dopamine = d*ones<vec>(csntc->m_N_CHA);

    inp = W_env2CRX *( 
        0.2*cos(12*M_PI*t/ti)+.3 +
        0.2*cos(24*M_PI*t/ti)+.3 ) /2.;  
 
    csntc->step(inp, 
            bg, 
            dopamine);
    csntc->store(t);
    activity = csntc->m_cortex.m_output;

    readout_pot += delta*(-readout_pot + W_readout*activity );
    readout = outfun(readout_pot);

    err = readout - target.row(testtimeline[tc%trials]).t();


    //
    // LEARN
    //
    
    if (learn and t >(ti*(1/2.)) and t<(ti*(3/4.)))
    {

        {
            static double eta = 1.4;
            static double beta = 0.00001;


            vec &x_prev = activity_prev;
            vec &y = readout;
            mat &w = W_readout;

            double dec = dot(x_prev,x_prev) + dot(inp_prev,inp_prev) + beta;
            vec g = (1-delta)*err_prev - err;

            w += (eta/delta)* g * (x_prev/dec).t();
        }

    }

    inp_prev = inp;
    activity_prev = activity;
    err_prev = err;

    return readout;
}
