import ROOT as R
import copy, sys, os, inspect 

# Usage from histFactory/plots/HHAnalysis/ : ./../../build/createHistoWithMultiDraw.exe -d ../../samples.json -o outputDir generatePlots.py 
scriptDir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.append(scriptDir)
from basePlotter import *
includes = ["/home/fynu/bfrancois/scratch/framework/oct2015/CMSSW_7_4_15/src/cp3_llbb/HHTools/histFactory_hh/readMVA.h"]

plots = []

objects = "nominal"
systematic = "nominal"
## Plots for ht and n vertex
basePlotter = BasePlotter()
plots.extend(basePlotter.generatePlots())

for plot in plots : 
    print plot
    #print plot["name"]
    print " "

