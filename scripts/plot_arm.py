#!/usr/bin/env python

from math import *
from pylab import *
import ik

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-s','--save',
        help="Print pictures",
        action="store_true", default=False)
args = parser.parse_args()

N_JOINTS = 3
N_CHANNELS = 3
N_CORTICAL = 300
N_SHAPES =3
LEARN_WINDOW = 144
SCALE = True
COMPLEX = True


def trial(teach, rout, LEARN_WINDOW, N_SHAPES, bg="#777777",fg="#000000",shapes=[0], width=3 ) :

    STIME = rout.shape[1] 
    TTIME = STIME/N_SHAPES
    
    rout = (teach.max() - teach.min())*((rout-.2)/.6) + teach.min()
    rng = zeros([LEARN_WINDOW,N_SHAPES])
    
    for x in xrange(N_SHAPES):
        rng[:,x] = arange(TTIME/2,TTIME/2+LEARN_WINDOW) + x*TTIME 
    rng = (rng.T.astype("int"))

    story = zeros([2,STIME])
    tstory = zeros([2,STIME])

    arm = ik.Arm(n=N_JOINTS, L=(ones(N_JOINTS)*(1/float(N_JOINTS-N_JOINTS/2.))))

    for k in shapes :
        for t in rng[k,:] :
     
            q = teach[:,t]
            xy = arm.get_xy(q) 
            tstory[:,t] = xy
        
            q = rout[:,t]
            xy = arm.get_xy(q) 
            story[:,t] = xy
        
        plot(*tstory[:,rng[k,:] ],markeredgewidth=4,
                    markeredgecolor=bg, 
                    linewidth=width,color=bg, alpha=.5 )
        plot(*story[:,rng[k,:] ],markeredgewidth=1,
                    markeredgecolor=fg, 
                    linewidth=1,color=fg )


def generalization() :
    
    data = loadtxt('out')
    rout  = data[(N_CORTICAL+6*N_CHANNELS+1*N_JOINTS):(N_CORTICAL+6*N_CHANNELS+2*N_JOINTS),:];
    
    teach1 = loadtxt("parameters/TEACH_1").T
    teach2 = loadtxt("parameters/TEACH_2").T
    teach3 = loadtxt("parameters/TEACH_3").T
    teach4 = loadtxt("parameters/TEACH_4").T

    teach=hstack([teach1,teach2,teach3,teach4])

    ion()

    STIME = N_SHAPES*480

    for shape in xrange(N_SHAPES):

        fig = figure(figsize=(2.5,2.5) )
        ax = fig.add_subplot(111)

        for x in xrange(N_SHAPES):
            trial(teach[:,(arange(STIME)+STIME*x)], 
                    rout[:,(arange(STIME)+STIME*x)],
                    LEARN_WINDOW, N_SHAPES,bg="#bbbbbb",fg="#666666",
                    shapes=[shape] )

        trial(teach[:,(arange(STIME)+STIME*(N_SHAPES))], 
                rout[:,(arange(STIME)+STIME*(N_SHAPES))],
                LEARN_WINDOW, N_SHAPES, shapes=[shape],width=4  )

        if SCALE :
            if COMPLEX :
                xlim([.6+.1*shape,1.4+.1*shape])
            else: 
                xlim([.8+.1*shape,1.2+.1*shape])
            ylim([.8,1.2])
        else :
            xlim([.55+.1*shape,1.45+.1*shape])
            ylim([.65,1.35])

        ax.set_aspect("equal")

def simpletest(shape=0) :
    
    data = loadtxt('out')
    rout  = data[(N_CORTICAL+6*N_CHANNELS+1*N_JOINTS):(N_CORTICAL+6*N_CHANNELS+2*N_JOINTS),:];
    teach = loadtxt("parameters/TEACH_1").T

    ion()
    
    fig = figure(figsize=(3.5,3.5))
    ax = fig.add_subplot(111)

    trial(teach,rout,LEARN_WINDOW, N_SHAPES,shapes=[shape])

    ax.set_aspect("equal")


# SIMPLE 4 SHAPES

# N_SHAPES=4
# N_CHANNELS=4
# simpletest(0)
# xlim([0.5,1.5])
# ylim([0.5,1.5])
# simpletest(1)
# xlim([0.5,1.5])
# ylim([0.5,1.5])
# simpletest(2)
# xlim([0.5,1.5])
# ylim([0.5,1.5])
# simpletest(3)
# xlim([0.5,1.5])
# ylim([0.5,1.5])
# raw_input()


generalization()
raw_input()
