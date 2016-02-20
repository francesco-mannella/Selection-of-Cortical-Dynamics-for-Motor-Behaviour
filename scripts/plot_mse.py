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

import sys
from pylab import *
ion()

#########################################################################
## MANAGE ARGUMENTS #####################################################
#########################################################################

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f','--file',
        help="file to read for data",
        action="store", default='mse_1')
parser.add_argument('-g','--gap',
        help="gap between series",
        action="store", default='20')
args = parser.parse_args()

filename = args.file
gap = int(args.gap)

while True :
    try :
        data = loadtxt(filename)
        if data.size != 0 :
            data = data[4:]
            fig = figure("mse")
            fig.clear()
            plot(data)
            xticks([0,data.size],array([0,data.size])*gap )
            fig.canvas.draw()
            pause(1)
    except :
        pass