# Selection-of-Cortical-Dynamics-for-Motor-Behaviour
This code implements the simulations described in the paper:
```
Selection of cortical dynamics for motor behaviour by the basal ganglia.
Mannella F, Baldassarre G.
Biol Cybern. 2015 Dec;109(6):575-95. doi: 10.1007/s00422-015-0662-6. Epub 2015 Nov 4.
PMID: 26537483 
```

You have open access to the paper at this [link](http://link.springer.com/article/10.1007%2Fs00422-015-0662-6).

## Simulations
The source code compiles three different simulations:
* Simulation 1 (simulation_loop_model_single_shape): a simulation of LOOP_MODEL learning to reproduce three different shapes, described in section 5.1 of the paper, and shown in the supplementary material 1.
* Simulation 2 (simulation_loop_model_multiple_shape): a simulation of LOOP_MODEL learning to generalize three different shapes over different scales or positions, also described in section 5.1 of the paper.
* Simulation 2 (simulation_loop_model_hand): a simulation of LOOP_MODEL learning to select and control three different positions in a 3D 20-DoF hand, described in section 5.2 of the paper and shown in  the supplementary material 2.

## Dependencies
This software depends on:
* The [Armadillo C++ linear algebra library](http://arma.sourceforge.net/docs.html)

It also requires
* The modules [regex](http://www.boost.org/doc/libs/1_60_0/libs/regex/doc/html/index.html), [system](http://www.boost.org/doc/libs/1_60_0/libs/system/doc/index.html), and [Filesystem](http://www.boost.org/doc/libs/1_60_0/libs/system/doc/index.html) from the [Boost C++ Libraries](http://www.boost.org/)

You finally need the [CENSLIB - Computational Embodied Neuroscience Simulation Library](http://censlib.sourceforge.net/) to compile simulation 3  

## Installation

* Download the zipped package from [here](https://github.com/francesco-mannella/Selection-of-Cortical-Dynamics-for-Motor-Behaviour/archive/master.zip)
* run these lines in the shell:
```
unzip Selection-of-Cortical-Dynamics-for-Motor-Behaviour.zip
cd Selection-of-Cortical-Dynamics-for-Motor-Behaviour
mkdir build && cd build
cmake .. && make -j 8 all
```
