#!/usr/bin/env python

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

from pylab import *
import matplotlib.gridspec as gridspec
import threading
import time

ion()
close('all')


def replot() :
    try :
        graphs = 6
        curr = 0

        data=loadtxt('out')
        mse=loadtxt('out_mse')

        data_cortex = data[:300,:];
        data_bg = data[-27:-9,:];
        data_rout = data[-6:-3,:];
        data_teach = data[-9:-6,:];
        data_rout = data[-6:-3,:];
        data_bginp = data[-3:,:];

        STIME=data.shape[1]

        gs = gridspec.GridSpec(3, 1)
        fig = figure("render",figsize=(14,8))
        fig.clear()

        subplot(gs[curr,0])
        plot(data_teach[0,:].T)
        plot(data_rout[0,:].T)
        ylim([0,1])
        xlim([1,STIME])

        curr += 1
        subplot(gs[curr,0])
        plot(data_teach[1,:].T)
        plot(data_rout[1,:].T)
        ylim([0,1])
        xlim([1,STIME])

        curr += 1
        subplot(gs[curr,0])
        plot(data_teach[2,:].T)
        plot(data_rout[2,:].T)
        ylim([0,1])
        xlim([1,STIME])

        fig.canvas.draw() 

        fig = figure("mse",figsize=(6,4))
        fig.clear()
        plot(mse)
        fig.canvas.draw() 
    

    except :
        pass

        
class KeyThread(threading.Thread):
    def run(self):
        global end_thread;
        
        raw_input()
        
        lock = threading.Lock()        
        lock.acquire()
        end_thread = True;
        lock.release()
   

end_thread = False
  
k = KeyThread()
k.start()

while(not end_thread) :
    replot()
    time.sleep(0.005)

k.join()