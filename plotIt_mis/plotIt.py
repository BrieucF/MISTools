#! /usr/bin/env python

import os, sys
import datetime, time

#usage python plotit.py [suffix]

plotNameAppend = str(datetime.date.today())#+"_"+str(int(round(time.time() * 1000)))
dirForPlots = "./"
suffix = ""
if len(sys.argv) == 2 : 
    suffix = "_"+sys.argv[1]

plotDir = dirForPlots+"plots_all_"+plotNameAppend+suffix
if not os.path.isdir(plotDir) :
    os.system("mkdir "+ plotDir)
print "../../plotIt/plotIt -o " + plotDir + " mis_plotter_all.yml"
os.system("../../plotIt/plotIt -o " + plotDir + " mis_plotter_all.yml -v -y")
