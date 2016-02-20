#!/usr/bin/env python

from math import *
from pylab import *
import sys
import os
import time
import ik

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-s','--save',
        help="Print pictures",
        action="store_true", default=False)
args = parser.parse_args()

import matplotlib.gridspec as gridspec
graphs = 3
gs = gridspec.GridSpec(5, 1)


MOON = 0
SQUARE = 1
INF = 2


data=loadtxt('out')
data_da = loadtxt('out_da')
data_inp = loadtxt('out_inp')

FSPYRAL = loadtxt("parameters/SPYRAL")
DDATA = [FSPYRAL,FSPYRAL,FSPYRAL]

labels = ['a','b','c']

TCH = hstack(DDATA)

TEACH_1 = loadtxt("parameters/TEACH_1").T
TEACH_2 = loadtxt("parameters/TEACH_2").T
TEACH_3 = loadtxt("parameters/TEACH_3").T
TEACH_4 = loadtxt("parameters/TEACH_4").T

TEACH = vstack([TEACH_1,TEACH_2,TEACH_3,TEACH_4])

teach = (TCH.max() - TCH.min())*((TEACH-.2)/.6) + TCH.min()

rout = data[-6:-3,:];
rout = (TCH.max() - TCH.min())*((rout-.2)/.6) + TCH.min()


#########################################################################
## PRESETS ##############################################################
#########################################################################

close('all')
ion()

timesteps=rout.shape[1]
timestep_gap = timesteps/3

# phase ranges
gp = 5/32.
rng1 = arange((timestep_gap*gp),(timestep_gap-timestep_gap*gp))
rng2 = arange((timestep_gap*gp),(timestep_gap-timestep_gap*gp)) + timestep_gap
rng3 = arange((timestep_gap*gp),(timestep_gap-timestep_gap*gp)) + 2*timestep_gap

trng1 = arange((timestep_gap/2),(timestep_gap-timestep_gap/4))
trng2 = arange((timestep_gap/2),(timestep_gap-timestep_gap/4)) + timestep_gap
trng3 = arange((timestep_gap/2),(timestep_gap-timestep_gap/4)) + 2*timestep_gap

rng = vstack([rng1,rng2,rng3])
trng = vstack([trng1,trng2,trng3])

n=3
arm = ik.Arm(n=n, L=ones(n)*(1/float(n-n/2.)))

xxgap = TCH.max()-TCH.min()

for CURR in xrange(3) :

    rng_curr = rng[CURR,:].astype("int") 

    tstory = zeros([rng_curr.size,2])
    story = zeros([rng_curr.size,2])
      
    
    for t in zip(trng[CURR,:]):
        
        q = TEACH[:,t]
        xy = arm.get_xy(q) 
        js = arm.get_joint_positions(q)
        tstory=roll(tstory.T,-1).T
        tstory[-1:,:] = xy
        tst = tstory[find(logical_and(tstory[:,0]!=0,tstory[:,1]!=0 )),:]

       
    for t,k in zip(rng_curr, range(rng_curr.size) ):
        
        q = rout[:,t]
        xy = arm.get_xy(q) 
        js = arm.get_joint_positions(q)
        story=roll(story.T,-1).T
        story[-1:,:] = xy
        st = story[find(logical_and(story[:,0]!=0,story[:,1]!=0 )),:]

        f = figure("post",figsize=(8,7))
        f.clear()
    
        plot(tst[:,0],tst[:,1], linewidth=.3,alpha=.4,color="b" ) 
        plot(st[-50::,0],st[-50::,1], linewidth=.8,color="r" ) 
        plot(st[-25::,0],st[-25::,1], linewidth=1.6,color="r" ) 
        plot(st[-10::,0],st[-10::,1], linewidth=2.6,color="r" ) 
        plot(st[-5::,0],st[-5::,1], linewidth=3.6,color="r" ) 
        for j in range(1,n+1) :
            plot( [ js[0,j-1], js[0,j] ], \
                    [ js[1,j-1], js[1,j] ], \
                    "o-",markeredgewidth=8,
                    markeredgecolor="#aaaaaa", 
                    linewidth=8,color="#aaaaaa" )
            plot( [ js[0,j-1], js[0,j] ], \
                    [ js[1,j-1], js[1,j] ], \
                    "o-",markeredgewidth=1,
                    markeredgecolor="#000000", 
                    linewidth=1, color="#000000"  ) 
        plot([ js[0,-1:], xy[0] ], \
                [ js[1,-1:], xy[1] ], \
                "o-",markeredgewidth=8,
                markeredgecolor="#aaaaaa", 
                linewidth=8,color="#aaaaaa" )
        plot([ js[0,-1:], xy[0] ], \
                [ js[1,-1:], xy[1] ], \
                "o-",markeredgewidth=2,
                markeredgecolor="#000000", 
                linewidth=1,color="#000000" ) 
        xlim([-.5,2])
        ylim([-.5,2])
        f.canvas.draw()
        if args.save == True :
            savefig("x_{:s}_{:05d}.png".format(labels[CURR],k))
     
        curr = 1
        prng = arange(t-50,t).astype("int") 
        
        fig = figure("render",figsize=(8,2))
        fig.clear()
        
        curr+=2
        subplot(gs[:curr,0])
        for x in xrange(300) :
            plot(data[x,:])
        xlim([prng[0],prng[-1]])
        xticks([])
        yticks([])
        axis("off")

        subplot(gs[curr,0])
        plot(data_da, color="#cccc33")
        fill_between(arange(timesteps),data_da,0,facecolor="#cccc33")
        ylim([-.2,1.2])
        xlim([prng[0],prng[-1]])
        xticks([])
        yticks([])
        axis("off")
        
        curr += 1
        subplot(gs[curr,0])
        for x in xrange(3) : 
            ccolor = zeros(3)
            ccolor[x]+=1
            fill_between(arange(timesteps),data[324+x,:],0,facecolor=ccolor)
        ylim([0,1])
        xlim([prng[0],prng[-1]])
        xticks([])
        yticks([])
        axis("off")
        fig.canvas.draw() 
        
        if args.save == True :
            savefig("x_acti_{:s}_{:05d}.png".format(labels[CURR],k),transparent=True )     

raw_input() 








