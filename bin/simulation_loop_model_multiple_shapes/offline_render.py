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

ion()
close('all')

import matplotlib.gridspec as gridspec

nj = 3
nb = 3
nc = 300

data=loadtxt('out')
data_da = loadtxt('out_da')
data_inp = loadtxt('out_inp')

data_cortex = data[:nc                         ,:];
data_bg = data[nc:(nc+6*nb     )               ,:];
data_teach = data[(nc+6*nb     ):(nc+6*nb+1*nj),:];
data_rout  = data[(nc+6*nb+1*nj):(nc+6*nb+2*nj),:];
data_bginp = data[(nc+6*nb+2*nj):              ,:];

print(data_rout.shape)
STIME=data.shape[1]




graphs = 6
curr = 1

gs = gridspec.GridSpec(14, 1)

fig = figure("render",figsize=(7,5))
fig.clear()

plot(data_rout.T)

 
curr+=4
subplot(gs[:curr,0])
ax=imshow(data_cortex,interpolation='none',aspect='auto',vmin=0,vmax=1)
cbaxes = fig.add_axes([0.92, 0.62, 0.01, 0.28]) 
colorbar(ax,cax=cbaxes)

curr += 2
subplot(gs[(curr-2):curr,0])
plot(data_teach[0,:].T)
plot(data_rout[0,:].T)
ylim([0,1])
xlim([1,STIME])

curr += 2
subplot(gs[(curr-2):curr,0])
plot(data_teach[1,:].T)
plot(data_rout[1,:].T)
ylim([0,1])
xlim([1,STIME])

curr += 2
subplot(gs[(curr-2):curr,0])
plot(data_teach[2,:].T)
plot(data_rout[2,:].T)
ylim([0,1])
xlim([1,STIME])

curr += 1
subplot(gs[curr,0])
plot(data_da)
ylim([-.2,1.2])
xlim([1,STIME])

curr += 1
subplot(gs[curr,0])
plot(data_bginp.T)
ylim([0,1])
xlim([1,STIME])
fig.canvas.draw() 


fig = figure("inp",figsize=(7,5))
plot(data_inp.T)
fig.canvas.draw() 

fig = figure("bg",figsize=(7,5))
imshow(data_bg, interpolation="none",aspect="auto")
fig.canvas.draw() 




raw_input()


