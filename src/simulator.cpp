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

const string Simulator::m_TEACH_FILE_1 = "parameters/TEACH_1";
const string Simulator::m_TEACH_FILE_2 = "parameters/TEACH_2";
const string Simulator::m_TEACH_FILE_3 = "parameters/TEACH_3";
const string Simulator::m_TEACH_FILE_4 = "parameters/TEACH_4";
const string Simulator::m_FUN_FILE_1 = "parameters/FUN_1";
const string Simulator::m_FUN_FILE_2 = "parameters/FUN_2";
const string Simulator::m_FUN_FILE_3 = "parameters/FUN_3";
const string Simulator::m_FUN_FILE_4 = "parameters/FUN_4";

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Simulator::Simulator(string name, int SEED): Parametrizable(name)
{

    // Use Parametrizable addParameter and saveParameters methods
    // to syncronize parameters

    addParameter("SESSIONS",m_SESSIONS);
    addParameter("FUN_PERIOD",m_FUN_PERIOD);
    addParameter("FUN_CENTER",m_FUN_CENTER);
    addParameter("FUN_PHASE",m_FUN_PHASE);
    addParameter("BG_BIAS_AMP",m_BG_BIAS_AMP);
    addParameter("BG_NOISE",m_BG_NOISE);
    addParameter("FUN_AMP",m_FUN_AMP);
    addParameter("DT",m_DT); 
    addParameter("TAU",m_TAU);
    addParameter("LAMBDA",m_LAMBDA);
    addParameter("BETA",m_BETA);
    addParameter("ETA",m_ETA);
    addParameter("LEARN_WINDOW",m_LEARN_WINDOW);
    addParameter("DA_GAP",m_DA_GAP);
    addParameter("PRINT",m_PRINT);
    addParameter("INP2BG",m_INP2BG);
    addParameter("INP2M",m_INP2M);

    try
    {
        loadParameters();
    }
    catch(parameter_file_exception &e)
    {

        m_SESSIONS = 2000;
        m_FUN_PERIOD = 30;
        m_FUN_CENTER = 0.5;
        m_FUN_PHASE = 0.5;
        m_FUN_AMP = 0.8;
        m_BG_BIAS_AMP = 0.5;
        m_BG_NOISE = 0.0;
        m_DT = 0.001; 
        m_TAU = 0.005;
        m_LAMBDA = 0.000001;
        m_ETA = 0.4;
        m_BETA = 0.00001;
        m_LEARN_WINDOW = 144;
        m_DA_GAP = 0.01;
        m_PRINT = "TRUE";
        m_INP2M = 1;
        m_INP2BG = 3; 
        saveParameters();

    }

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // Init m_csntc main object
    m_csntc = unique_ptr<CSNTC>(new CSNTC(name+"_csntc",SEED));   
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////   

    m_SELECTIONS = m_csntc->m_N_CHA;

    m_DELTA = m_DT/m_TAU;    // Init m_DELTA
    m_STIME = m_csntc->m_STIME;    // Init m_STIME
    m_TTIME = m_STIME/m_SELECTIONS;    //  M_TTIME is 1/3 of m_STIME


    // Define a sequence of the timestep indices of the whole session
    uvec m_arng = conv_to<uvec>::from(linspace(0,m_STIME-1,m_STIME));

    // Init learning windows' ranges
    //
    m_rng = zeros<uvec>(m_SELECTIONS*m_LEARN_WINDOW);
    m_rngs = zeros<umat>(m_LEARN_WINDOW,m_SELECTIONS);
    for(int x =0;  x<m_SELECTIONS; x++)
    {
        uvec rng = m_arng( span(m_TTIME/2,m_TTIME/2 +m_LEARN_WINDOW - 1) )+ x*m_TTIME; 
        m_rng( span(x*m_LEARN_WINDOW, (x+1)*m_LEARN_WINDOW - 1) ) = rng;
        m_rngs(span::all, x) = rng;
    }
    
    // Init trial sequence
    m_trial_idx = conv_to<uvec>::from(linspace(0,m_SELECTIONS-1,m_SELECTIONS));

    // Init cortical indices
    m_crx_idx = conv_to<uvec>::from(
            linspace(0,m_csntc->m_N_CRX-1,m_csntc->m_N_CRX)); 


    // Init input-to-cortex weights
    m_W_inp2m = m_INP2M*randu<mat>(m_csntc->m_N_CRX, m_csntc->m_N_CRX);
    m_W_inp2m( span(0, m_csntc->m_N_CRX*(1-m_csntc->m_THA2CRX_SPREAD) ), span::all) *= 0;

    // Init input-to-bg weights
    m_W_inp2bg = m_INP2BG;

    // Init dopamine
    //       ____________          _____________          _____________
    // |____/::::::::::::\___||___/:::::::::::::\___||___/:::::::::::::\___|
    // |                     ||                     ||                     |
    // |______1-st_trial_____||______2-nd trial_____||______3-rd trial_____|

    m_DA_GAP = ceil(m_DA_GAP/m_DT); 
    m_da = mod(linspace(0,m_STIME-1,m_STIME), m_TTIME);     
    m_da = conv_to<vec>::from( (m_da>m_DA_GAP)%(m_da<(m_TTIME-3*m_DA_GAP))%(m_da>(3*m_DA_GAP)) );     
    double min_val = 0.2;    // Minimum of dopamine 
    double increment = 0.6;  // Dopaminergic increment when switched on
    m_da = m_da*increment +min_val;  

    // Init cortical input and target output
    build_funcs(Simulator::m_TEACH_FILE_1,Simulator::m_FUN_FILE_1, /* reset */ true );

    modify_trial_sequence(true);  


}

/////////////////////////////////////////////////////////////////////////

void Simulator::build_funcs( const string &teach_sec_file, const string &fun_sec_file, bool reset_all )
{
    // CORTICAL INPUT                   |                      |                      |
    {
        mat fn;
        fn.load(fun_sec_file,raw_ascii);      
        
        if( remind(m_csntc->m_N_CRX, fn.n_rows) != 0 )
            throw(runtime_error("FUN dimension is not multiple of cortical N"));
        
        int N = m_csntc->m_N_CRX/fn.n_rows;
        m_fun = repmat(fn,N,1);
    }

    // TARGET OUTPUT
    mat teaches;
    teaches.load(teach_sec_file,raw_ascii);    // Take data from file

    double min_val = 0.2;    // Minimum 
    double scale = 0.6;     // Range to maximum

    // Normalize between target between (min_val) and (min_val+range)
    teaches = min_val +scale*( (teaches -teaches.min()) / (teaches.max()-teaches.min() ) );

    if (reset_all)
    {
        // The number of readout units is the number of columns from the target output data
        m_n_readout = teaches.n_cols;

        // Init vectors readout and target vectors 
        m_rout = zeros<mat>(m_n_readout, m_STIME); 
        m_Wout = randn<mat>(m_n_readout, m_csntc->m_N_CRX)/m_csntc->m_N_CRX;     
        m_teach = zeros<mat>(m_n_readout, m_STIME);
    }

    // Fill m_teach with the value of target output only within the learning windows
    for(int teach=0; teach<m_n_readout; teach++)
    {
        uvec uteach = uidx(teach);
        for(int x=0; x<m_SELECTIONS; x++)
        
            m_teach(uteach,m_rngs.col(x)) = teaches( m_rngs.col(x), uteach ).t();     
    }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Simulator::modify_trial_sequence(bool reset)
{

    // Shuffle the trial sequence
    if(reset==true)
        m_trial_idx = conv_to<uvec>::from(linspace(0,m_SELECTIONS-1,m_SELECTIONS));
    else
        m_trial_idx = shuffle(m_trial_idx);

    ////////////////////////////////////////////////////////

    // Each unit of the reservoir gets the input signal
    m_crx_inp = m_fun;

    ////////////////////////////////////////////////////////

    // Each striatum unit gets a noise input
    mat bgm = randn<mat>(m_SELECTIONS,m_STIME)*m_BG_NOISE;

    for(int s=0; s<m_SELECTIONS; s++ )
    {
        // Only one striatal unit is biased during each trial 
        bgm(m_trial_idx(s),span(s*m_TTIME,(s+1)*m_TTIME-1)) += m_BG_BIAS_AMP;
    }
    // Negative values are truncated
    bgm = bgm%(bgm>0);

    // Updating the final matrix
    m_bg_inp = zeros<mat>(m_SELECTIONS,m_STIME);
    m_bg_inp +=  bgm; 

    ////////////////////////////////////////////////////////


    // Shuffling the trials of target outputs
    m_curr_teach = zeros<mat>(m_n_readout, m_STIME);
    for(int teach=0; teach< m_n_readout; teach++)
    {
        uvec uteach = uidx(teach);
        for(int s=0; s<m_SELECTIONS; s++ )
            m_curr_teach( uteach, m_rngs.col(s) ) = m_teach( uteach , m_rngs.col( m_trial_idx(s)) );
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Simulator::step(int t)
{
    // A simulation step
    m_csntc->step(  
            m_W_inp2m*m_crx_inp.col(t),     // Input to cortex
            m_W_inp2bg*m_bg_inp.col(t),     // Input to striatum
            ones<vec>(m_csntc->m_N_CHA)*m_da(t) );     // Dopamine

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Simulator::run_offline()
{
    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////
    // Cleaning

    // Pre-training sessions
    for( int session=0; session<1; session++ )
    {
        // Shuffle trial on each session 
        modify_trial_sequence(true);  

        // Session
        for( int t=0; t<m_STIME; t++ )
        {
            step(t);   // Spread
        }
    }

    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////
    // Learning readout

    // Pre-learning sessions
    for( int session=0; session<3; session++ )
    {

        // Shuffle trial on each session 
        modify_trial_sequence(true);

        for( int t=0; t<m_STIME; t++ )
        {

            step(t);  // Spread    
            m_csntc->store(t);  // Store

            // Spread of readout units
            {
                vec &x = m_csntc->m_cortex.m_output;
                mat &W = m_Wout;
                for(int out=0; out< m_n_readout; out++)
                    m_rout(out,t) = dot( W.row(out), x );
            }

        }
    }

    // Batch ridge regressions
    {
        mat X = m_csntc->m_cortex.m_data.cols(m_rng);
        int N = m_csntc->m_cortex.m_N;

        for(int out=0; out< m_n_readout; out++)
        {
            uvec uout = uidx(out);
            vec Y = m_curr_teach(uout, m_rng).t();
            m_Wout.row(out) = (inv( X * X.t() + m_LAMBDA*eye(N,N) )* (X*Y)).t() ; 
        }
    }

    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////
    // Test 

    // Used to store striatum inputs during the test session
    mat data_inps = zeros<mat>(m_csntc->m_N_CHA,m_STIME); 

    // Shuffle trials
    modify_trial_sequence(true);
    // Test session
    for( int t=0; t<m_STIME; t++ )
    {
        step(t);    // Spread 
        m_csntc->store(t);    // Store

        // Store inputs to the striatum
        data_inps(span::all,t) = (m_W_inp2bg*m_bg_inp.col(t));

        // Spread of readout units
        {
            vec &x = m_csntc->m_cortex.m_output;
            mat &W = m_Wout;
            for(int out=0; out< m_n_readout; out++)
                m_rout(out,t) = dot( W.row(out), x );
        }
    }

    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////
    // Saving

    {
        if ( m_PRINT == "TRUE" )
        {
            // Save main data in "out"
            mat dataout;

            dataout =                     m_csntc->m_cortex.m_data;
            dataout = join_cols(dataout,  m_csntc->m_thalamus.m_data);
            dataout = join_cols(dataout,  m_csntc->m_ganglia.m_data);
            dataout = join_cols(dataout,  m_curr_teach );
            dataout = join_cols(dataout,  m_rout );
            dataout = join_cols(dataout,  data_inps);

            stringstream sout(""); sout << "out"; 
            ofstream fout("out"); 
            dataout.raw_print(fout);

            // Save da activity in "out_da"
            ofstream fda("out_da");
            m_da.raw_print(fda);

            // Save cortical input in "out_inp"
            ofstream finp("out_inp");
            m_fun.raw_print(finp);

        }

        // Test is selections were optimal during the test session 
        bool comp = test_bg_competition(
                m_csntc->m_thalamus.m_data,
                m_rngs);

        // Compute the mean of NRMSE in the three trials of the sessions 
        double mse_tot = 0;
        for(int out=0; out< m_n_readout; out++)
        {
            uvec uout = uidx(out);
            mse_tot +=  mse(m_curr_teach(uout,m_rng) , m_rout(uout,m_rng))/m_n_readout;
        }

        // Store mean NRMSE and seletion test in "mse"
        ofstream fmse("mse");

        fmse << mse_tot  <<  " " << comp <<  endl;

    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Simulator::run_offline_double_teach()
{
    // TRAINING
    {
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        // TEACH1

        build_funcs(m_TEACH_FILE_1,m_FUN_FILE_1);

        ///////////////////////////////////////////////////////////
        // Cleaning

        // Pre-training sessions
        for( int session=0; session<1; session++ )
        {
            // Shuffle trial on each session 
            modify_trial_sequence(true);  

            // Session
            for( int t=0; t<m_STIME; t++ )
            {
                step(t);   // Spread
            }
        }

        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        // Learning readout

        // Pre-learning sessions
        for( int session=0; session<3; session++ )
        {

            // Shuffle trial on each session 
            modify_trial_sequence(true);  

            for( int t=0; t<m_STIME; t++ )
            {

                step(t);  // Spread    
                m_csntc->store(t);  // Store

                // Spread of readout units
                {
                    vec &x = m_csntc->m_cortex.m_output;
                    mat &W = m_Wout;
                    for(int out=0; out< m_n_readout; out++)
                        m_rout(out,t) = dot( W.row(out), x );
                }

            }
        }

        mat cortex_all_data_1 = m_csntc->m_cortex.m_data;
        mat curr_teach_all_1 = m_curr_teach;

        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        //TEACH_2

        build_funcs(m_TEACH_FILE_2,m_FUN_FILE_2);

        ///////////////////////////////////////////////////////////
        // Cleaning

        // Pre-training sessions
        for( int session=0; session<1; session++ )
        {
            // Shuffle trial on each session 
            modify_trial_sequence(true);  

            // Session
            for( int t=0; t<m_STIME; t++ )
            {
                step(t);   // Spread
            }
        }

        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        // Learning readout 

        // Pre-learning sessions
        for( int session=0; session<3; session++ )
        {

            // Shuffle trial on each session 
            modify_trial_sequence(true);  

            for( int t=0; t<m_STIME; t++ )
            {

                step(t);  // Spread    
                m_csntc->store(t);  // Store

                // Spread of readout units
                {
                    vec &x = m_csntc->m_cortex.m_output;
                    mat &W = m_Wout;
                    for(int out=0; out< m_n_readout; out++)
                        m_rout(out,t) = dot( W.row(out), x );
                }

            }
        }

        mat cortex_all_data_2 = m_csntc->m_cortex.m_data;
        mat curr_teach_all_2 = m_curr_teach;
        
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        //TEACH_3

        build_funcs(m_TEACH_FILE_3,m_FUN_FILE_3);

        ///////////////////////////////////////////////////////////
        // Cleaning

        // Pre-training sessions
        for( int session=0; session<1; session++ )
        {
            // Shuffle trial on each session 
            modify_trial_sequence(true);  

            // Session
            for( int t=0; t<m_STIME; t++ )
            {
                step(t);   // Spread
            }
        }

        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        // Learning readout 

        // Pre-learning sessions
        for( int session=0; session<3; session++ )
        {

            // Shuffle trial on each session 
            modify_trial_sequence(true);  

            for( int t=0; t<m_STIME; t++ )
            {

                step(t);  // Spread    
                m_csntc->store(t);  // Store

                // Spread of readout units
                {
                    vec &x = m_csntc->m_cortex.m_output;
                    mat &W = m_Wout;
                    for(int out=0; out< m_n_readout; out++)
                        m_rout(out,t) = dot( W.row(out), x );
                }

            }
        }

        mat cortex_all_data_3 = m_csntc->m_cortex.m_data;
        mat curr_teach_all_3 = m_curr_teach;


        mat cortex_all_data = join_rows(
                join_rows(cortex_all_data_1,cortex_all_data_2),cortex_all_data_3);
        mat curr_teach_all = join_rows(
                join_rows(curr_teach_all_1,curr_teach_all_2),curr_teach_all_3);


        // Batch ridge regressions
        {
            mat X = join_rows( 
                    join_rows( cortex_all_data_1.cols(m_rng),
                        cortex_all_data_2.cols(m_rng)),
                    cortex_all_data_3.cols(m_rng));       

            int N = m_csntc->m_cortex.m_N;

            for(int out=0; out< m_n_readout; out++)
            {
                uvec uout = uidx(out);
                vec Y = join_rows(
                        join_rows( curr_teach_all_1(uout, m_rng),
                            curr_teach_all_2(uout, m_rng)), 
                        curr_teach_all_3(uout, m_rng)  ) .t();

                m_Wout.row(out) = ( inv( X * X.t() + m_LAMBDA*eye(N,N) )* (X*Y) ).t() ; 
            }
        }
    } 

    // TEST
    {
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        // TEST1 

        // Used to store striatum inputs during the test session
        mat data_inps = zeros<mat>(m_csntc->m_N_CHA,m_STIME); 

        ///////////////////////////////////////////////////////////   

        build_funcs(m_TEACH_FILE_1,m_FUN_FILE_1);

        ///////////////////////////////////////////////////////////
        // Cleaning
        // Pre-training sessions
        for( int session=0; session<1; session++ )
        {
            // Shuffle trial on each session 
            modify_trial_sequence(true);  

            // Session
            for( int t=0; t<m_STIME; t++ )
            {
                step(t);   // Spread
            }
        }
        ///////////////////////////////////////////////////////////

        // Shuffle trials
        modify_trial_sequence(true);  
        
        // Test session
        for( int t=0; t<m_STIME; t++ )
        {
            step(t);    // Spread 
            m_csntc->store(t);    // Store

            // Store inputs to the striatum
            data_inps(span::all,t) = (m_W_inp2bg*m_bg_inp.col(t));

            // Spread of readout units
            {
                vec &x = m_csntc->m_cortex.m_output;
                mat &W = m_Wout;
                for(int out=0; out< m_n_readout; out++)
                    m_rout(out,t) = dot( W.row(out), x );
            }
        }

        mat rout_1 = m_rout; 
        mat data_inps_1 = data_inps; 

        mat cortex_all_data_1 = m_csntc->m_cortex.m_data;
        mat thalamus_all_data_1 = m_csntc->m_thalamus.m_data;
        mat ganglia_all_data_1 = m_csntc->m_ganglia.m_data;
        mat curr_teach_all_1 = m_curr_teach;
        mat fun_all_1 = m_fun;

        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        // TEST2


        build_funcs(m_TEACH_FILE_2,m_FUN_FILE_2);

        ///////////////////////////////////////////////////////////
        // Cleaning
        // Pre-training sessions
        for( int session=0; session<1; session++ )
        {
            // Shuffle trial on each session 
            modify_trial_sequence(true);  

            // Session
            for( int t=0; t<m_STIME; t++ )
            {
                step(t);   // Spread
            }
        }
        ///////////////////////////////////////////////////////////


        // Shuffle trials
        modify_trial_sequence(true);  

        // Test session
        for( int t=0; t<m_STIME; t++ )
        {
            step(t);    // Spread 
            m_csntc->store(t);    // Store

            // Store inputs to the striatum
            data_inps(span::all,t) = (m_W_inp2bg*m_bg_inp.col(t));

            // Spread of readout units
            {
                vec &x = m_csntc->m_cortex.m_output;
                mat &W = m_Wout;
                for(int out=0; out< m_n_readout; out++)
                    m_rout(out,t) = dot( W.row(out), x );
            }
        }



        mat rout_2 = m_rout; 
        mat data_inps_2 = data_inps; 

        mat cortex_all_data_2 = m_csntc->m_cortex.m_data;
        mat thalamus_all_data_2 = m_csntc->m_thalamus.m_data;
        mat ganglia_all_data_2 = m_csntc->m_ganglia.m_data;
        mat curr_teach_all_2 = m_curr_teach;
        mat fun_all_2 = m_fun;

        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        // TEST3 


        build_funcs(m_TEACH_FILE_3,m_FUN_FILE_3);

        ///////////////////////////////////////////////////////////
        // Cleaning
        // Pre-training sessions
        for( int session=0; session<1; session++ )
        {
            // Shuffle trial on each session 
            modify_trial_sequence(true);  

            // Session
            for( int t=0; t<m_STIME; t++ )
            {
                step(t);   // Spread
            }
        }
        ///////////////////////////////////////////////////////////


        // Shuffle trials
        modify_trial_sequence(true);  

        // Test session
        for( int t=0; t<m_STIME; t++ )
        {
            step(t);    // Spread 
            m_csntc->store(t);    // Store

            // Store inputs to the striatum
            data_inps(span::all,t) = (m_W_inp2bg*m_bg_inp.col(t));

            // Spread of readout units
            {
                vec &x = m_csntc->m_cortex.m_output;
                mat &W = m_Wout;
                for(int out=0; out< m_n_readout; out++)
                    m_rout(out,t) = dot( W.row(out), x );
            }
        }

   
        mat rout_3 = m_rout; 
        mat data_inps_3 = data_inps; 

        mat cortex_all_data_3 = m_csntc->m_cortex.m_data;
        mat thalamus_all_data_3 = m_csntc->m_thalamus.m_data;
        mat ganglia_all_data_3 = m_csntc->m_ganglia.m_data;
        mat curr_teach_all_3 = m_curr_teach;
        mat fun_all_3 = m_fun;
       
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        // TEST4 


        build_funcs(m_TEACH_FILE_4,m_FUN_FILE_4);

        ///////////////////////////////////////////////////////////
        // Cleaning
        // Pre-training sessions
        for( int session=0; session<1; session++ )
        {
            // Shuffle trial on each session 
            modify_trial_sequence(true);  

            // Session
            for( int t=0; t<m_STIME; t++ )
            {
                step(t);   // Spread
            }
        }
        ///////////////////////////////////////////////////////////


        // Shuffle trials
        modify_trial_sequence(true);  

        // Test session
        for( int t=0; t<m_STIME; t++ )
        {
            step(t);    // Spread 
            m_csntc->store(t);    // Store

            // Store inputs to the striatum
            data_inps(span::all,t) = (m_W_inp2bg*m_bg_inp.col(t));

            // Spread of readout units
            {
                vec &x = m_csntc->m_cortex.m_output;
                mat &W = m_Wout;
                for(int out=0; out< m_n_readout; out++)
                    m_rout(out,t) = dot( W.row(out), x );
            }
        }

   
        mat rout_4 = m_rout; 
        mat data_inps_4 = data_inps; 

        mat cortex_all_data_4 = m_csntc->m_cortex.m_data;
        mat thalamus_all_data_4 = m_csntc->m_thalamus.m_data;
        mat ganglia_all_data_4 = m_csntc->m_ganglia.m_data;
        mat curr_teach_all_4 = m_curr_teach;
        mat fun_all_4 = m_fun;


        ///////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////
        // Saving
     
        mat cortex_all_data = 
            join_rows(
                join_rows(
                    join_rows(cortex_all_data_1,
                        cortex_all_data_2),
                    cortex_all_data_3),
                cortex_all_data_4);
        mat thalamus_all_data = 
            join_rows(
                join_rows(
                    join_rows(thalamus_all_data_1,
                        thalamus_all_data_2),
                    thalamus_all_data_3),
                thalamus_all_data_4);
        mat ganglia_all_data = 
            join_rows(
                join_rows(
                    join_rows(ganglia_all_data_1,
                        ganglia_all_data_2),
                    ganglia_all_data_3),
                ganglia_all_data_4);
        mat curr_teach_all = 
            join_rows(
                join_rows(
                    join_rows(curr_teach_all_1,
                        curr_teach_all_2),
                    curr_teach_all_3),
                curr_teach_all_4);
        mat fun_all = 
            join_rows(
                join_rows(
                    join_rows(fun_all_1,
                        fun_all_2),
                    fun_all_3),
                fun_all_4);



        mat rout_all = 
            join_rows(
                    join_rows(
                        join_rows(rout_1,
                            rout_2),
                        rout_3),
                    rout_4);
        mat inps_all = 
            join_rows(
                    join_rows(
                        join_rows(data_inps_1,
                            data_inps_2),
                        data_inps_3),
                    data_inps_4);

        {
            if ( m_PRINT == "TRUE" )
            {
                // Save main data in "out"
                mat dataout;

                dataout =                     cortex_all_data;
                dataout = join_cols(dataout,  thalamus_all_data);
                dataout = join_cols(dataout,  ganglia_all_data);
                dataout = join_cols(dataout,  curr_teach_all );
                dataout = join_cols(dataout,  rout_all );
                dataout = join_cols(dataout,  inps_all);

                stringstream sout(""); sout << "out"; 
                ofstream fout("out"); 
                dataout.raw_print(fout);

                // Save da activity in "out_da"
                ofstream fda("out_da");
                join_cols(join_cols(m_da,m_da),m_da).raw_print(fda);

                // Save cortical input in "out_inp"
                ofstream finp("out_inp");
                fun_all.raw_print(finp);

            }

            // Test is selections were optimal during the test session 
            bool comp = test_bg_competition(
                    m_csntc->m_thalamus.m_data,
                    m_rngs);

            // Compute the mean of NRMSE in the three trials of the sessions 
            double mse_tot = 0;
            for(int out=0; out< m_n_readout; out++)
            {
                uvec uout = uidx(out);
                mse_tot +=  mse(m_curr_teach(uout,m_rng) , m_rout(uout,m_rng))/m_n_readout;
            }

            // Store mean NRMSE and seletion test in "mse"
            ofstream fmse("mse");
            fmse << mse_tot  <<  " " << comp <<  endl;

        }
    }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Simulator::run_online()
{
    // Prepare online variables
    Layer readout_layer("csntc_readout");
    vec err_prev = zeros<vec>(m_n_readout);
    vec x_prev = zeros<vec>(m_csntc->m_N_CRX);
    vec inp_prev = zeros<vec>(m_csntc->m_N_CRX);

    { ofstream fmse("out_mse");  }

    // Pre-learning sessions
    for( int session=0; session<m_SESSIONS; session++ )
    {

        // Used to store striatum inputs during the test session
        mat data_inps = zeros<mat>(m_csntc->m_N_CHA,m_STIME); 

        // Shuffle trial on each session 
        modify_trial_sequence();

        bool is_test_session = (session % 20 == 0); 

        for( int t=0; t<m_STIME; t++ )
        {

            step(t);  // Spread    

            {

                // Prepare local variables
                vec &x = m_csntc->m_cortex.m_output;
                mat &W = m_Wout;  
                readout_layer.step(W*x);    // Readout spreading

                if (not is_test_session)   // LEARNING
                {
                    // Prepare local variables
                    vec &y =  readout_layer.m_output;
                    vec ty =  m_curr_teach( span::all, t );
                    vec err = y - ty; 

                    // Prepare decorrelation and derivative
                    int DRVTV = 0;     // Switch off computing the derivative of error instead 
                    vec d = x_prev / ( dot(inp_prev,inp_prev) + dot(x_prev,x_prev) + m_BETA ); 
                    vec g = DRVTV*(1 - m_DELTA)*err_prev - err;

                    // Update 
                    W +=  (sum(ty)>0)*(m_ETA/m_DELTA)*( g*d.t() );

                    // Update stores
                    err_prev = err;
                    x_prev = x;
                    inp_prev = m_crx_inp.col(t);

                }
                else if(is_test_session)    // TEST
                {
                    // Print on test 
                    m_csntc->store(t);   
                    readout_layer.store(t);

                    // Store inputs to the striatum
                    data_inps(span::all,t) = (m_W_inp2bg*m_bg_inp.col(t));
                }

            }
        }

        if (is_test_session)
        {

            if ( m_PRINT == "TRUE" )
            {
                // Save main data in "out"
                mat dataout;

                dataout =                     m_csntc->m_cortex.m_data;
                dataout = join_cols(dataout,  m_csntc->m_thalamus.m_data);
                dataout = join_cols(dataout,  m_csntc->m_ganglia.m_data);
                dataout = join_cols(dataout,  m_curr_teach );
                dataout = join_cols(dataout,  readout_layer.m_data );
                dataout = join_cols(dataout,  data_inps);

                stringstream sout(""); sout << "out"; 
                ofstream fout("out"); 
                dataout.raw_print(fout);

                // Save da activity in "out_da"
                ofstream fda("out_da");
                m_da.raw_print(fda);

                // Save cortical input in "out_inp"
                ofstream finp("out_inp");
                m_fun.raw_print(finp);

            }

            // Compute the mean of NRMSE in the three trials of the sessions 
            double mse_tot = 0;
            for(int out=0; out< m_n_readout; out++)
            {
                uvec uout = uidx(out);
                mse_tot +=  mse(m_curr_teach(uout,m_rng) , readout_layer.m_data(uout,m_rng))/m_n_readout;
            }
            ofstream fmse("out_mse",fstream::app);
            fmse << mse_tot << endl;

        }
    }
}


