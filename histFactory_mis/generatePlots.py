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

basePlotter = BasePlotter()
weightEpression = "((runOnMC)? (btaggingReweightingMM*leptonsReweightingweight*lumiReweightingLumiWeight*mcReweightingweight*((m_dataset.cross_section==3503.71)? ((abs(jetmetbjet1Flavor)==5 && abs(jetmetbjet2Flavor)==5 && jetmetnj==2)*1.148+(abs(jetmetbjet1Flavor)==5 && abs(jetmetbjet2Flavor)==5 && jetmetnj>2)*1.348+ ( (abs(jetmetbjet1Flavor)!=5 && abs(jetmetbjet2Flavor)==5) || (abs(jetmetbjet1Flavor)==5 && abs(jetmetbjet2Flavor)!=5) )*1.348+( abs(jetmetbjet1Flavor)!=5 && abs(jetmetbjet2Flavor)!=5 )*1.359) : 1.)) : 1.)"
plots.extend(basePlotter.generatePlots(categories = ["All"], stage = "no_cut", requested_plots = ["bdtoutput"], weightExpression = weightEpression, extraCut = "", systematic = "nominal", extraString = ""))

for plot in plots : 
    print plot
    #print plot["name"]
    print " "

