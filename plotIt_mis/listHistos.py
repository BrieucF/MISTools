#! /usr/bin/env python

import yaml

#usage : python listHisto.py [yields]
from ROOT import * 
import argparse

parser = argparse.ArgumentParser(description='Facility to produce the yml with plots information.')
parser.add_argument('--yields', help='If you just want to produce the yields and systematics.', action="store_true")
parser.add_argument('--dyamc', help='If you just want to produce the Mll histograms.', action="store_true")
parser.add_argument('-m', '--mass', dest='mass', default="650", help='Mass point for the yields.')
parser.add_argument('-d', '--directory', dest='directory', default="combinedBDT_v0", help='Directory of the input rootfiles.')
parser.add_argument('--blinded', dest='unblinded', help='If you want to produce blinded plots', action="store_false")
parser.add_argument('--mjj', help='Use if you want to produce plots related to Mjj', action="store_true")
args = parser.parse_args()

baseDir = "/home/fynu/bfrancois/scratch/framework/oct2015/CMSSW_7_4_15/src/cp3_llbb/" 
fileName = baseDir+"/MISTools/histFactory_mis/"+args.directory+"/TT_histos.root"

skim = False

file = TFile(fileName) 
keys = file.GetListOfKeys() 
alreadyIn = []

# Dictionary containing all the plots
plots = {}

logY = 'both'
for_yield = False
if args.yields :
    logY = False
    for_yield = True
    print "Will produce yields"
defaultStyle = {
        'log-y': logY,
        'save-extensions': ['pdf', 'png'],
        'legend-columns': 2,
        'show-ratio': True,
        'show-overflow': True,
        'show-errors': True,
        'for-yields': for_yield
        }

defaultStyle_events_per_gev = defaultStyle.copy()
defaultStyle_events_per_gev.update({
        'y-axis': 'Events',
        'y-axis-format': '%1% / %2$.2f GeV',
        })

defaultStyle_events = defaultStyle.copy()
defaultStyle_events.update({
        'y-axis': 'Events',
        'y-axis-format': '%1% / %2$.2f',
        })
print args

for key in keys :
    if key.GetName() not in alreadyIn  and not "__" in key.GetName():

        if not "jj_M_" in key.GetName() and args.mjj : 
            continue

        alreadyIn.append(key.GetName())
        plot = {
                'x-axis': key.GetName(),
                }

        if "lep1_pt" in key.GetName() :
            plot['x-axis'] = "Leading lepton p_{T} (GeV)"
            plot.update(defaultStyle_events_per_gev)
        elif "lep2_pt" in key.GetName() :  
            plot['x-axis'] = "Sub-leading lepton p_{T} (GeV)"
            plot.update(defaultStyle_events_per_gev)
        elif "jet1_pt" in key.GetName() :  
            plot['x-axis'] = "Leading jet p_{T} (GeV)"
            plot.update(defaultStyle_events_per_gev)
        elif "jet2_pt" in key.GetName() :  
            plot['x-axis'] = "Sub-leading jet p_{T} (GeV)"
            plot.update(defaultStyle_events_per_gev)
        elif "lep1_eta" in key.GetName() :  
            plot['x-axis'] = "Leading jet #eta"
            plot.update(defaultStyle_events)
        elif "lep2_eta" in key.GetName() :  
            plot['x-axis'] = "Sub-leading jet #eta"
            plot.update(defaultStyle_events)
        elif "jet1_eta" in key.GetName() :  
            plot['x-axis'] = "Leading jet #eta"
            plot.update(defaultStyle_events)
        elif "jet2_eta" in key.GetName() :  
            plot['x-axis'] = "Sub-leading jet #eta"
            plot.update(defaultStyle_events)
        elif "lep1_phi" in key.GetName() :  
            plot['x-axis'] = "Leading lepton #phi"
            plot.update(defaultStyle_events)
        elif "lep2_phi" in key.GetName() :  
            plot['x-axis'] = "Sub-leading lepton #phi"
            plot.update(defaultStyle_events)
        elif "jet1_phi" in key.GetName() :  
            plot['x-axis'] = "Leading jet #phi"
            plot.update(defaultStyle_events)
        elif "jet2_phi" in key.GetName() :  
            plot['x-axis'] = "Sub-leading jet #phi"
            plot.update(defaultStyle_events)
        elif "jet1_CSV" in key.GetName() :  
            plot['x-axis'] = "Leading jet CSVv2 discriminant"
            plot.update(defaultStyle_events)
        elif "jet2_CSV" in key.GetName() :  
            plot['x-axis'] = "Sub-leading jet CSVv2 discriminant"
            plot.update(defaultStyle_events)
        elif "jet1_JP" in key.GetName() :  
            plot['x-axis'] = "Leading jet JP discriminant"
            plot.update(defaultStyle_events)
        elif "jet2_JP" in key.GetName() :  
            plot['x-axis'] = "Sub-leading jet JP discriminant"
            plot.update(defaultStyle_events)
        elif "ll_M_" in key.GetName() :  
            plot['x-axis'] = "m_{ll} (GeV)"
            plot.update(defaultStyle_events_per_gev)
        elif "ll_pt_" in key.GetName() :  
            plot['x-axis'] = "Dilepton system p_{T} (GeV)"
            plot.update(defaultStyle_events_per_gev)
        elif "jj_pt_" in key.GetName() :  
            plot['x-axis'] = "Dijet system p_{T} (GeV)"
            plot.update(defaultStyle_events_per_gev)
        elif "met_pt" in key.GetName() :  
            plot['x-axis'] = "#slash{E}_{T} (GeV)"
            plot.update(defaultStyle_events_per_gev)
        elif "met_phi" in key.GetName() :  
            plot['x-axis'] = "#phi_{#slash{E}_{T}}"
            plot.update(defaultStyle_events)
        elif "ll_DR_l_l_All_hh_llmetjj_HWWleptons_btagM_csv_cleaning_cut" in key.GetName() :  
            plot['x-axis'] = "#DeltaR(leading lepton, sub-leading lepton)"
            plot.update(defaultStyle_events)
        elif "ll_DR_l_l" in key.GetName() :  
            plot['x-axis'] = "#DeltaR(leading lepton, sub-leading lepton)"
            plot.update(defaultStyle_events)
        elif "jj_DR_j_j" in key.GetName() :  
            plot['x-axis'] = "#DeltaR(leading jet, sub-leading jet)"
            plot.update(defaultStyle_events)
        elif "ll_DPhi_l_l" in key.GetName() :  
            plot['x-axis'] = "#Delta#phi(leading lepton, sub-leading lepton)"
            plot.update(defaultStyle_events)
        elif "jj_DPhi_j_j" in key.GetName() :  
            plot['x-axis'] = "#Delta#phi(leading jet, sub-leading jet)"
            plot.update(defaultStyle_events)
        elif "llmetjj_pt_" in key.GetName() :  
            plot['x-axis'] = "p_{T}^{lljj#slash{E}_{T}}"
            plot.update(defaultStyle_events_per_gev)
        elif "llmetjj_M_" in key.GetName() :  
            plot['x-axis'] = "m_{lljj#slash{E}_{T}}"
            plot.update(defaultStyle_events_per_gev)
        elif "DPhi_ll_met_" in key.GetName() :  
            plot['x-axis'] = "#Delta#phi(ll, #slash{E}_{T})"
            plot.update(defaultStyle_events)
        elif "DPhi_ll_jj" in key.GetName() :  
            plot['x-axis'] = "#Delta#phi(ll, jj)"
            plot.update(defaultStyle_events)
        elif "minDPhi_l_met_" in key.GetName() :  
            plot['x-axis'] = "min(#Delta#phi(lepton, #slash{E}_{T}))"
            plot.update(defaultStyle_events)
        elif "maxDPhi_l_met_" in key.GetName() :  
            plot['x-axis'] = "max(#Delta#phi(lepton, #slash{E}_{T}))"
            plot.update(defaultStyle_events)
        elif "MT_" in key.GetName() :  
            plot['x-axis'] = "m_{ll#slash{E}_{T}}"
            plot.update(defaultStyle_events_per_gev)
        elif "MTformula_" in key.GetName() :  
            plot['x-axis'] = "MT"
            plot.update(defaultStyle_events_per_gev)
        elif "projMET_" in key.GetName() :  
            plot['x-axis'] = "Projected #slash{E}_{T}"
            plot.update(defaultStyle_events_per_gev)
        elif "DPhi_jj_met" in key.GetName() :  
            plot['x-axis'] = "#Delta#phi(jj, #slash{E}_{T})"
            plot.update(defaultStyle_events)
        elif "minDPhi_j_met" in key.GetName() :  
            plot['x-axis'] = "min#Delta#phi(j, #slash{E}_{T})"
            plot.update(defaultStyle_events)
        elif "maxDPhi_j_met" in key.GetName() :  
            plot['x-axis'] = "max#Delta#phi(j, #slash{E}_{T})"
            plot.update(defaultStyle_events)
        elif "minDR_l_j" in key.GetName() :  
            plot['x-axis'] = "min#DeltaR(l, j)"
            plot.update(defaultStyle_events)
        elif "maxDR_l_j" in key.GetName() :  
            plot['x-axis'] = "max#DeltaR(l, j)"
            plot.update(defaultStyle_events)
        elif "DR_ll_jj_" in key.GetName() :  
            plot['x-axis'] = "#DeltaR(ll, jj)"
            plot.update(defaultStyle_events)
        elif "DR_llmet_jj" in key.GetName() :  
            plot['x-axis'] = "#DeltaR(ll#slash{E}_{T}, jj)"
            plot.update(defaultStyle_events)
        elif "DPhi_ll_jj_" in key.GetName() :  
            plot['x-axis'] = "#DeltaPhi(ll, jj)"
            plot.update(defaultStyle_events)
        elif "nllmetjj_" in key.GetName() :  
            plot['x-axis'] = "#llmetjj"
            plot.update(defaultStyle_events)
        elif "nLep_" in key.GetName() :  
            plot['x-axis'] = "Number of leptons"
            plot.update(defaultStyle_events)
        elif "nJet_" in key.GetName() :  
            plot['x-axis'] = "Number of jets"
            plot.update(defaultStyle_events)
        elif "nBJetMediumCSV_" in key.GetName() :  
            plot['x-axis'] = "Number of b-tagged jets (CSVv2 medium)"
            plot.update(defaultStyle_events)
        elif "cosThetaStar_CS" in key.GetName() :  
            plot['x-axis'] = "cos(#theta^{*}_{CS})"
            plot.update(defaultStyle_events)
        elif "MVA" in key.GetName() and "mjj_cr" in key.GetName() :
            plot['x-axis'] = "BDT output"
            plot.update(defaultStyle_events)
            if not args.unblinded :
                plot['blinded-range'] = [0, 0.6]
        elif "jj_M_" in key.GetName() and "_cr_" in key.GetName() :  
            plot['x-axis'] = "m_{jj} (GeV)"
            plot.update(defaultStyle_events_per_gev)
            if not args.unblinded :
                plot['blinded-range'] = [75, 140]
        elif "MVA_" in key.GetName() :
            plot['x-axis'] = "BDT output"
            plot.update(defaultStyle_events)
            if not args.unblinded :
                plot['blinded-range'] = [0, 0.6]
                plot['for-yields'] = True
        elif "jj_M_" in key.GetName() :  
            plot['x-axis'] = "m_{jj} (GeV)"
            plot.update(defaultStyle_events_per_gev)
            if not args.unblinded :
                plot['blinded-range'] = [75, 140]
        elif "isElEl" in key.GetName() and "highBDT_mjjP_" + args.mass in key.GetName() :
            plot['x-axis'] = "isElEl"
            plot.update(defaultStyle_events)
            plot['for-yields'] = True
            plot['yields-title'] = "high-BDT %s, m$_{jj}$-P"%(args.mass)
            plot['yields-table-order'] = 0
        elif "isElEl" in key.GetName() and "highBDT_mjjSB_" + args.mass in key.GetName() :
            plot['x-axis'] = "isElEl"
            plot.update(defaultStyle_events)
            plot['for-yields'] = True
            plot['yields-title'] = "high-BDT %s, m$_{jj}$-SB"%(args.mass)
            plot['yields-table-order'] = 1
        elif "isElEl" in key.GetName() and "lowBDT_mjjP_" + args.mass in key.GetName() :
            plot['x-axis'] = "isElEl"
            plot.update(defaultStyle_events)
            plot['for-yields'] = True
            plot['yields-title'] = "low-BDT %s, m$_{jj}$-P"%(args.mass)
            plot['yields-table-order'] = 2
        elif "isElEl" in key.GetName() and "lowBDT_mjjSB_" + args.mass in key.GetName() :
            plot['x-axis'] = "isElEl"
            plot.update(defaultStyle_events)
            plot['for-yields'] = True
            plot['yields-title'] = "low-BDT %s, m$_{jj}$-SB"%(args.mass)
            plot['yields-table-order'] = 3
        else:
            plot.update(defaultStyle_events)

        if "gen_" in key.GetName():
            flavour = key.GetName().split("_")[1]
            plot['x-axis'] = "is "+flavour
            plot['no-data'] = True
            plot.update(defaultStyle_events)
        if any(x in key.GetName() for x in ['sr_400_ext', 'sr_650_ext']) and not args.unblinded :
            plot['no-data'] = True
        if "scaleFactor" in key.GetName() :
            plot['x-axis'] = "Scale factor"
            plot.update(defaultStyle_events)
            plot['no-data'] = True

        label_x = 0.22
        label_y = 0.895
        if "MuMu" in key.GetName() : 
            plot['labels'] = [{
                'text': '#mu#mu channel',
                'position': [label_x, label_y],
                'size': 24
                }]
        elif "MuEl" in key.GetName() : 
            plot['labels'] = [{
                'text': '#mue + e#mu channels',
                'position': [label_x, label_y],
                'size': 24
                }]
        elif "All" in key.GetName() : 
            plot['labels'] = [{
                'text': '#mu#mu + ee + #mue + e#mu channels',
                'position': [label_x, label_y],
                'size': 24
                }]
        elif "ElEl" in key.GetName() : #need to be the last one as we use isElEl to compute the yields 
            plot['labels'] = [{
                'text': 'ee channel',
                'position': [label_x, label_y],
                'size': 24
                }]
        if args.yields and "isElEl" in key.GetName() and args.mass in key.GetName() :
            plots[key.GetName()] = plot
        if args.dyamc and "ll_M_" in key.GetName() :
            plots[key.GetName()] = plot
        if not args.yields and not args.dyamc :
            plots[key.GetName()] = plot

with open("allPlots.yml", "w") as f:
    yaml.dump(plots, f)
