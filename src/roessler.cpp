#!/usr/bin/python 

# CSNTC - Cortico-striato-nigro-thalamo-cortical comutational model
# Copyright (C) 2014 Francesco Mannella <francesco.mannella@gmail.com>
#
# This file is part of CSNTC.
#
# CSNTC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CSNTC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CSNTC.  If not, see <http://www.gnu.org/licenses/>.

from math import *
from pylab import *

from scipy.ndimage import gaussian_filter1d
import scipy.optimize  
import time 

close('all')

class Arm:
     
    def __init__(self, n=3,q=None, q0=None, L=None):
        """Set up the basic parameters of the arm.
        All lists are in order [shoulder, elbow, wrist].
         
        :param list q: the initial joint angles of the arm
        :param list q0: the default (resting state) joint configuration
        :param list L: the arm segment lengths
        """
        # initial joint angles
        if q is None: q = zeros(n)
        self.q = q
        # some default arm positions
        if q0 is None: q0 = ones(n)* pi*2
        self.q0 = q0
        # arm segment lengths
        if L is None: L = ones(n-1) 
        self.L = L
        
        self.n = n
        self.max_angles = ones(n)* pi
        self.min_angles = ones(n)* -pi/8.
 
    def get_xy(self, q=None):
        if q is None: q = self.q
        x = sum([ self.L[j]*cos( q[0:(j+1)].sum() ) for j in range(q.size) ])
        y = sum([ self.L[j]*sin( q[0:(j+1)].sum() ) for j in range(q.size) ])
        return [x, y]
 
    def get_joint_positions(self,q):
        """This method finds the (x,y) coordinates of each joint"""
        
        x = array([  sum([ self.L[j]*cos( q[0:(k+1)].sum() ) for k in range(j) ]) for j in range(self.n) ])
        x = hstack([0,x])

 
        y = array([ sum( [ self.L[j]*sin(q[0:(k+1)].sum() ) for k in range(j) ]) for j in range(self.n) ])
        y = hstack([0,y])


        return array([x, y])
    
    def inv_kin(self, xy):
 
        def distance_to_default(q, *args): 
            # weights found with trial and error, get some wrist bend, but not much
            weight = ones(self.n)
            weight = array([0.1,0.2,.2])
            return sqrt(sum([(qi - q0i)**2 * wi for qi,q0i,wi in zip(q, self.q0, weight)]))
 
        def x_constraint(q, xy):
            x = sum([ self.L[j]*cos( q[0:(j+1)].sum() ) for j in range(q.size) ]) -xy[0]
            return x
 
        def y_constraint(q, xy): 
            y = sum([ self.L[j]*sin( q[0:(j+1)].sum() ) for j in range(q.size) ]) -xy[1]
            return y
        
        return scipy.optimize.fmin_slsqp( func=distance_to_default, 
            x0=self.q, eqcons=[x_constraint, y_constraint], 
            args= (xy,), bounds= array([self.min_angles,self.max_angles]).T, iprint=0) # iprint=0 suppresses output

             

if __name__ == "__main__" : 
    """A function for plotting an arm, and having it calculate the 
    inverse kinematics such that given the position it 
    finds the appropriate joint angles to reach that point."""
 
    MOON = 0
    SQUARE = 1
    INF = 2
    CIRCLE = 3
    TRIANGLE = 4
    SPYRAL = 5
    NAME = ""

    PLOT_H = True
 
    ion()

    n = 3
    T = 48
    TN = 10

    # TRANSLATION
    #                 -----MOON--------    -----SQUARE----    -----INF-----      
    params = array([ [ 0.8,  1.0,   0.2,   0.8,   1.0, 0.2,   0.8, 1.0, 0.2 ],  # TEACH_1
                     [ 0.2,  1.5,   0.2,   1.5,   0.2, 0.2,   1.2, 0.2, 0.2 ],  # TEACH_2
                     [ 0.5,  1.25,  0.2,   1.15,  0.6, 0.2,   1.0, 0.6, 0.2 ],  # TEACH_3
                     [ 0.35, 1.375, 0.2,   1.325, 0.4, 0.2,   1.1, 0.4, 0.2 ]])  # TEACH_4
    
    # scale
    #                 -----MOON--------    -----SQUARE----    -----INF-----      
    params = array([ [ 0.8,  1.0,   0.05,   0.8,   1.0, 0.05,   0.8, 1.0, 0.05 ],  # TEACH_1
                     [ 0.8,  1.0,   0.1,    0.8,   1.0, 0.1,    0.8, 1.0, 0.1  ],  # TEACH_2
                     [ 0.8,  1.0,   0.2,    0.8,   1.0, 0.2,    0.8, 1.0, 0.2  ],  # TEACH_3
                     [ 0.8,  1.0,   0.15,   0.8,   1.0, 0.15,   0.8, 1.0, 0.15 ]])  # TEACH_4


    shapes = {'MOON':0,'SQUARE':1,'INF':2 } 
    teachs = {'TEACH_1':0,'TEACH_2':1,'TEACH_3':2, 'TEACH_4':3 } 
    
    Nshapes = len(shapes.keys())

    # create an instance of the arm
    arm = Arm(n=n, L=ones(n)*(1/float(n-n/2.)))
    
    for teach in ['TEACH_1','TEACH_2','TEACH_3', 'TEACH_4'] :
        current_teach = zeros([T*TN*Nshapes, n])
        for shape in ['MOON','SQUARE','INF'] :
            b =   params[ teachs[teach], 3*shapes[shape] ]  
            c =   params[ teachs[teach], 3*shapes[shape] + 1 ]  
            r =   params[ teachs[teach], 3*shapes[shape] + 2 ]  

            X  = array([])
            Y  = array([])
            
            CURRENT = shapes[shape] 


            if CURRENT == SPYRAL :
            
                NAME="SPYRAL"
                
                rr = linspace(r,0,T)

                XX = b + rr*cos(linspace(0,6*pi,T))
                YY = c + rr*sin(linspace(0,6*pi,T))
                X = XX 
                Y = YY 

            if CURRENT == CIRCLE :
            
                NAME="CIRCLE"

                XX = b + r*cos(linspace(0,2*pi,T))
                YY = c + r*sin(linspace(0,2*pi,T))
                X = XX 
                Y = YY  

            if CURRENT == INF :
            
                NAME="INF"

                XX0 = b + r*cos(linspace(0,2*pi,T/2))
                YY0 = c + r*sin(linspace(0,2*pi,T/2))
                XX1 = b +2*r+ r*cos(linspace(pi-pi/float(T/2),-pi,T/2))
                YY1 = c +r*sin(linspace(pi-pi/float(T/2),-pi,T/2))
                X = hstack([XX0,XX1])
                Y = hstack([YY0,YY1])

            if CURRENT == SQUARE :

                NAME="SQUARE"

                # square
                ##### points
                P =zeros([4,2])
                P[0,:] = [b + r*cos(0), c + r*sin(0) ] 
                P[1,:] = [b + r*cos(pi/2), c + r*sin(pi/2) ] 
                P[2,:] = [b + r*cos(pi), c + r*sin(pi) ] 
                P[3,:] = [b + r*cos(pi*(3/2.)), c + r*sin(pi*(3/2.)) ] 

                
                for k in range(1) :
                    for t in range(4) :
                        
                        tt = t%4 
                        ttt = (t+1)%4 
                        
                        mmax = array([P[tt,0],P[ttt,0]]).max()
                        mmin = array([P[tt,0],P[ttt,0]]).min()
                        XX = linspace( mmin,mmax ,T/4 )
                        
                        if P[tt,0].max() == mmax :
                            XX = XX[::-1] 
                        
                        YY = \
                                ( ( P[ttt,1]-P[tt,1] )/( P[ttt,0]-P[tt,0] )  )*\
                                (XX - P[tt,0]) + P[tt,1]
                        X = hstack([X,XX])
                        Y = hstack([Y,YY])
            
            if CURRENT == MOON :

                NAME="MOON"
                XX0 = c + r*cos(linspace(-pi*2/3.,pi*2/3.,T/2))
                YY0 = b + r*sin(linspace(-pi*2/3.,pi*2/3.,T/2))
                
                r1 = r*sin(pi*2/3.)
                c1 = c+r*cos(pi*2/3.)
                XX1 = c1 + r1*cos(linspace(pi/2.,-pi/2.,T/2))
                YY1 = b + r1*sin(linspace(pi/2.,-pi/2.,T/2))
            
            
                X = hstack([XX0,XX1])
                Y = hstack([YY0,YY1])


            if CURRENT == TRIANGLE :
            
                NAME="TRIANGLE"

                # triangle
                ##### points
                P =zeros([3,2])
                P[0,:] = [b + r*cos(-pi/4.), c + r*sin(-pi/4.) ] 
                P[1,:] = [b + r*cos(pi/2), c + r*sin(pi/2) ] 
                P[2,:] = [b + r*cos(pi+pi/4.), c + r*sin(pi+pi/4.) ] 
            
                for k in range(1) :
                    for t in range(3) :
                        
                        tt = t%3 
                        ttt = (t+1)%3 
                        
                        mmax = array([P[tt,0],P[ttt,0]]).max()
                        mmin = array([P[tt,0],P[ttt,0]]).min()
                        XX = linspace( mmin,mmax ,T/3 )
                        
                        if P[tt,0].max() == mmax :
                            XX = XX[::-1] 
                        
                        YY = \
                                ( ( P[ttt,1]-P[tt,1] )/( P[ttt,0]-P[tt,0] )  )*\
                                (XX - P[tt,0]) + P[tt,1]
                        X = hstack([X,XX])
                        Y = hstack([Y,YY])
        
            T = X.size
            
            print T
            
            #close("all")

            dataxy = zeros([2,T])
            dataq = zeros([n,T])
            datajxy = zeros([2,n+1,T])

        
            for t in range(T):
                pos = array([X[t], Y[t]])
                q = arm.inv_kin(xy=pos) # get new arm angles
                xy = arm.get_xy(q)
                error = np.sqrt(( array([X[t],Y[t]]) - array(xy))**2)
                if norm(error) < 1e-6 :
                    dataxy[:,t] = xy
                    dataq[:,t] = q
                    
                    js = arm.get_joint_positions(q)
                    datajxy[:,:,t] = js

                    if PLOT_H :
                        f = figure("pre",figsize=(8,7))
                        f.clear()
                        plot(X, Y, linewidth=4 ) 
                        plot( X[t], Y[t], scalex = 1, scaley = 1)
                        for j in range(1,n+1) :
                            plot([ js[0,j-1], js[0,j] ], \
                                    [ js[1,j-1], js[1,j] ], \
                                    "o-",markeredgewidth=8, linewidth=8 )
                        plot([ js[0,-1:], xy[0] ], \
                                [ js[1,-1:], xy[1] ], \
                                "o-",markeredgewidth=8, linewidth=8 )
                        plot( xy[0], xy[1], "o")
                        xlim([-.5,2])
                        ylim([-.5,2])
                        f.canvas.draw()
                else :
                    dataxy[:,t] = None
                    dataq[:,t] = None
                    datajxy[:,:,t] = None

            
            # interpolate NaNs
            for qq in dataq :
                ok = -np.isnan(qq)
                xp = ok.ravel().nonzero()[0]
                fp = qq[-np.isnan(qq)]
                x  = isnan(qq).ravel().nonzero()[0]
                qq[np.isnan(qq)] = interp(x, xp, fp)
            for qq in dataxy :
                ok = -np.isnan(qq)
                xp = ok.ravel().nonzero()[0]
                fp = qq[-np.isnan(qq)]
                x  = isnan(qq).ravel().nonzero()[0]
                qq[np.isnan(qq)] = interp(x, xp, fp)
        

            f = figure("curves",figsize=(16,4))
            subplot(131)
            plot(dataxy.T)
            subplot(132)
            plot(dataq.T)
            subplot(133)
            plot(dataxy[0],dataxy[1])
            xlim([-.5,2])
            ylim([-.5,2])
            
            if PLOT_H :
                for t in range(T):
                    f = figure("post",figsize=(8,7))
                    pos = array([X[t], Y[t]])
                    q = dataq[:,t]
                    xy = arm.get_xy(q)
                    error = np.sqrt(( array([X[t],Y[t]]) - array(xy))**2)
                    
                    js = arm.get_joint_positions(q)

                    f.clear()
                    plot(X, Y, linewidth=4 ) 
                    plot( X[t], Y[t], scalex = 1, scaley = 1)
                    for j in range(1,n+1) :
                        plot( [ js[0,j-1], js[0,j] ], \
                                [ js[1,j-1], js[1,j] ], \
                                "o-",markeredgewidth=8, linewidth=8 )
                    plot([ js[0,-1:], xy[0] ], \
                            [ js[1,-1:], xy[1] ], \
                            "o-",markeredgewidth=8, linewidth=8 )
                    plot( xy[0], xy[1], "o")
                    xlim([-.5,2])
                    ylim([-.5,2])
                    f.canvas.draw()
                    

            dataq = repeat(dataq.T,TN,0)

            start = (T*TN*shapes[shape])
            current_teach[ start:(start + T*TN),:] = dataq

        savetxt(teach, current_teach)

             
