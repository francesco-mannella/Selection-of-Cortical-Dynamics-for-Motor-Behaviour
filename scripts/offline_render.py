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
#imshow(data_inp[298:,:], interpolation="none",aspect="auto")
fig.canvas.draw() 

fig = figure("bg",figsize=(7,5))
imshow(data_bg, interpolation="none",aspect="auto")
fig.canvas.draw() 




raw_input()

# 
# class bg:
#     N = 3
#     l_sd1 = arange(3) 
#     l_sd2 = arange(3)  + 3
#     l_stn = arange(3)  + 3*2
#     l_gpi = arange(3)  + 3*3
#     l_gpe = arange(3)  + 3*4
# 
# ################################################################################################
# ################################################################################################
# ################################################################################################
# ## Plot
# colors=['r','g','b']
# scale=.22
# rowscale=.27
# start=.0
# gap=0
# 
# stime=data.shape[1]
# inp=data[324:,:]
# 
# ylims=[-.1,1.1]
# 
# row =rowscale*tile(arange(bg.N)+1,(stime,1)).T 
# rng = arange(stime)[(stime*(1/3.)-stime*(1/9.)):(stime*(1/3.))+stime*(1/9.)] 
# 
# fig = figure(facecolor='white',figsize=(3.8,5))
# subplots_adjust(
#     left=0.01, right=.99, 
#     top=.99, bottom=.01,
#     wspace=.001,hspace=.1) 
# 
# d = inp 
# p=subplot(9,1,1)
# for q in range(bg.N) :
#     qrng = rng -q*gap
#     plot( qrng ,d[q,rng] ,lw=2,c=colors[q],alpha=.5)
#     fill_between(qrng,d[q,rng],0,facecolor=colors[q], alpha=.1 )
# ylim(ylims)
# 
# p.set_axis_off()
# 
# d = data_da
# p=subplot(9,1,2)
# plot( qrng ,d[rng] ,lw=2,c="#cccc22",alpha=.5)
# fill_between(qrng,d[rng],0,facecolor="#cccc22", alpha=.1 )
# ylim(ylims)
# p.set_axis_off()
# 
# for (num,label) in enumerate([bg.l_sd1,bg.l_sd2,bg.l_stn,bg.l_gpi,bg.l_gpe]) :
#     p=subplot(9,1,num+3)
#     d = data_bg[label,:] 
#     for q in range(bg.N) :
#         qrng = rng -q*gap
#         plot( qrng ,d[q,rng] ,lw=2,c=colors[q],alpha=.5)
#         fill_between(qrng,d[q,rng],0,facecolor=colors[q], alpha=.1 )
#     ylim(ylims)
#     p.set_axis_off()
# 
# d = data_bg[:3,:]
# p=subplot(9,1,8)
# for q in range(bg.N) :
#     qrng = rng -q*gap
#     plot( qrng ,d[q,rng] ,lw=2,c=colors[q],alpha=.5)
#     fill_between(qrng,d[q,rng],0,facecolor=colors[q], alpha=.1 )
# 
# ylim(ylims)
# p.set_axis_off()
# 
# p=subplot(9,1,9)
# 
# 
# for x in xrange(10) :
#     d = data[(arange(3)*50)+25 + int(rand()*10)-5 ,:]
#     for q in range(bg.N) :
#         qrng = rng -q*gap
#         plot( qrng ,d[q,rng] ,lw=2,c=colors[q],alpha=.1)
#         fill_between(qrng,d[q,rng],0,facecolor=colors[q], alpha=.05 )
# 
# ylim(ylims)
# 
# 
# p.set_axis_off()
# 
# fig.canvas.draw()
# 
# 
# raw_input()
