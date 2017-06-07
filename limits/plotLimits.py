#!/usr/bin/env python
import sys
import numpy as np
import json
import collections
# User imports
tmpargv = sys.argv[:] # to prevent pyroot to hijack argparse we need to go around
sys.argv = []
# ROOT imports
import ROOT
from ROOT import TChain, TGraph, TCanvas, TLatex, TLegend, TLine, TColor
ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch()
ROOT.gROOT.ProcessLine(".x setTDRStyle.C")
ROOT.TGaxis.SetMaxDigits(3)
ROOT.PyConfig.IgnoreCommandLineOptions = True
sys.argv = tmpargv
# Get on with pyroot
c1 = TCanvas()
color = TColor()

limitMethod = 'Asymptotic'
#limitMethod = 'HybridNew'
outfileSuffix = '%s_withoutSystematics' % limitMethod
extra = ''
if limitMethod == 'MarkovChainMC': # or limitMethod == 'HybridNew':
    extra = '.123456'
unblind = False
spins = [0, 2]
MVAspin = 0
fbUnits = True
scaledToHH = False
#pathname = "/home/fynu/obondu/Higgs/CMSSW_7_1_5/src/HHStatAnalysis/"
pathname = '/home/fynu/bfrancois/scratch/framework/combine/CMSSW_7_4_7/src/HHStatAnalysis'
lowmjjcut = 95
higmjjcut = 135
mvacut = 0.1
masses = [
    [260, 270, 300, 350, 400, 450], # MVAmass = 400
    [450, 500, 550, 600, 650, 700, 800, 900] # MVAmass = 650
    ]
MVAmasses = [400, 650]
transitionLines = [450]
DIR = 'CARDS'
#DIR = 'STUDIES_spin'
# NB: only this line is needed to switch on and off the interpolation
interpolation = True
mgranularity = 5
allmasses = [np.arange(min(masslist), max(masslist) + mgranularity, mgranularity) for masslist in masses]
allmasses = [masslist.tolist() for masslist in allmasses]
dumpForCommonPlot = False

theoryCurves = [
    [ # spin 0 curves
#        'theoryPredictions/radion_Lambda-0.4.json',
#        'theoryPredictions/radion_Lambda-1.0.json',
#        'theoryPredictions/2HDM_type-1_tanb-0.10_cosba-0.10_m12-0_mh-125.json',
#        'theoryPredictions/2HDM_type-1_tanb-0.10_cosba-0.05_m12-0_mh-125.json',
#        'theoryPredictions/2HDM_type-1_tanb-0.05_cosba-0.10_m12-0_mh-125.json',
#        'theoryPredictions/2HDM_type-1_tanb-0.05_cosba-0.05_m12-0_mh-125.json',
#        'theoryPredictions/2HDM_type-2_tanb-0.10_cosba-0.10_m12-0_mh-125.json',
#        'theoryPredictions/2HDM_type-2_tanb-0.10_cosba-0.05_m12-0_mh-125.json',
#        'theoryPredictions/2HDM_type-2_tanb-0.05_cosba-0.10_m12-0_mh-125.json',
#        'theoryPredictions/2HDM_type-2_tanb-0.05_cosba-0.05_m12-0_mh-125.json',
    ],
    [ # spin 2 curves
#        'theoryPredictions/RS_KKGraviton_ktilda-0.3_kL-35.json',
    ],
]

BRHTo2B = 5.824e-1
BRHToVVTo2L2Nu = 1.055e-2
BRbbVV = 2 * BRHTo2B * BRHToVVTo2L2Nu


masses_ = masses
if interpolation:
    masses_ = allmasses

for ispin, spin in enumerate(spins):
    exp_2sigLow = {}
    exp_1sigLow = {}
    exp = {}
    exp_1sigHig = {}
    exp_2sigHig = {}
    obs = {}
    # TGraphs
    graph_obs = {}
    graph_obs_fullsim = {}
    graph_exp = {}
    graph_exp1sig = {}
    graph_exp2sig = {}
    # Various options
    c1.SetLogy(1)
    ymax = 1e4
    ymin = 1e1
    yaxisextra = ""
    if scaledToHH: 
        ymax = 1e6
        ymin = 1e3
    xmin = 200
    xmax = 950
    strategies = []
    for MVAmass in MVAmasses:
        strategies.append("cut-n-count_%i" % MVAmass)
    
    latex = True
    outfilename = 'limits_%s_spin-%i' % (outfileSuffix, spin)
    if scaledToHH:
        outfilename += '_scaledToHH'
    if latex:
        print("Latex dump in tables/%s.tex" % outfilename)
        with open('tables/%s.tex' % outfilename, 'w') as f:
            f.write(R'\begin{tabular}{@{}ccccc@{}} \toprule' + '\n')
            f.write(R'm$_X^{\mathrm{spin-}%i}$ (GeV) & Observed (fb) & Expected (fb) & 1 Standard deviation (fb) & 2 Standard deviations (fb) \\ \midrule' % spin + '\n')
    for iMVAmass, MVAmass in enumerate(MVAmasses):
        dumpForCommonPlotOutfile = 'tables/CMS_2b2w_low_mass_spin-%i.txt' % spin
        if MVAmass == 650:
            dumpForCommonPlotOutfile = 'tables/CMS_2b2w_high_mass_spin-%i.txt' % spin
        if dumpForCommonPlot:
            with open(dumpForCommonPlotOutfile, 'w') as dumpfile:
                dumpfile.write('mX\tobserved\texpected\n')
        for ikey, key in enumerate(strategies):
            if str(MVAmass) not in key:
                continue
            exp_2sigLow[key] = []
            exp_1sigLow[key] = []
            exp[key] = []
            exp_1sigHig[key] = []
            exp_2sigHig[key] = []
            obs[key] = []
        # READ THE FILES
            basename = "Test"
            if "cut-n-count" in key:
                basename = "CutBased"
            for mass in masses_[iMVAmass]:
                filename = basename + '_X%iToHHTo2B2L2Nu_BDT_X%i_%s_VS_TT_DY_MJJ-%s-%s_MVA-%s_All' % (spin, MVAspin, MVAmass, str(lowmjjcut), str(higmjjcut), str(mvacut).replace('.','p').replace('-','m'))
                chain = TChain("limit")
                chain.Add("%s%s/%i/higgsCombine%s.%s.mH%i%s.root" % (pathname, DIR, mass, filename, limitMethod, mass, extra))
                if limitMethod == 'Asymptotic':
                    for ievt in xrange(chain.GetEntries()):
                        chain.GetEntry( ievt )
                        l = chain.limit # default is in fb for sigma * BR
                        if not fbUnits:
                            l = l / 1.0e3
                        if scaledToHH:
                            l = l / BRbbVV
                        if ievt == 0:
                            exp_2sigLow[key].append( l )
                        elif ievt == 1:
                            exp_1sigLow[key].append( l )
                        elif ievt == 2:
                            exp[key].append( l )
                        elif ievt == 3:
                            exp_1sigHig[key].append( l )
                        elif ievt == 4:
                            exp_2sigHig[key].append( l )
                        elif ievt == 5:
                            obs[key].append( l )
                elif limitMethod == 'MarkovChainMC':
                    allexp = []
                    for ievt in xrange(chain.GetEntries()):
                        chain.GetEntry( ievt )
                        l = chain.limit # default is in fb for sigma * BR
                        if not fbUnits:
                            l = l / 1.0e3
                        if scaledToHH:
                            l = l / BRbbVV
                        allexp.append(l)
                    # NOTE: the following is not fully consistent with the output of combine which takes value from the toys intead of interpolation between toys for the percentile
                    allexp_np = allexp
                    exp_2sigHig[key].append(np.percentile(allexp_np, 97.5)) #, interpolation='higher') 
                    exp_1sigHig[key].append(np.percentile(allexp_np, 84.0)) #, interpolation='higher') 
                    exp[key].append(np.percentile(allexp_np, 50.0)) 
                    exp_1sigLow[key].append(np.percentile(allexp_np, 16.0)) #, interpolation='lower') 
                    exp_2sigLow[key].append(np.percentile(allexp_np,  2.5)) #, interpolation='lower') 
                elif limitMethod == 'HybridNew':
                    allexp = []
                    for ievt in xrange(chain.GetEntries()):
                        chain.GetEntry( ievt )
                        l = chain.limit # default is in fb for sigma * BR
                        if not fbUnits:
                            l = l / 1.0e3
                        if scaledToHH:
                            l = l / BRbbVV
                        allexp.append(l)
                    # NOTE: the following is not fully consistent with the output of combine which takes value from the toys intead of interpolation between toys for the percentile
                    allexp_np = allexp
                    exp_2sigHig[key].append(np.percentile(allexp_np, 97.5)) #, interpolation='higher') 
                    exp_1sigHig[key].append(np.percentile(allexp_np, 84.0)) #, interpolation='higher') 
                    exp[key].append(np.percentile(allexp_np, 50.0)) 
                    exp_1sigLow[key].append(np.percentile(allexp_np, 16.0)) #, interpolation='lower') 
                    exp_2sigLow[key].append(np.percentile(allexp_np,  2.5)) #, interpolation='lower') 
            
            graph_obs[key] = TGraph()
            graph_exp[key] = TGraph()
            graph_exp1sig[key] = TGraph()
            graph_exp2sig[key] = TGraph()
            exp_2sigLow[key] = np.asarray(exp_2sigLow[key], dtype='f')
            exp_1sigLow[key] = np.asarray(exp_1sigLow[key], dtype='f')
            exp[key] = np.asarray(exp[key], dtype='f')
            obs[key] = np.asarray(obs[key], dtype='f')
            exp_1sigHig[key] = np.asarray(exp_1sigHig[key], dtype='f')
            exp_2sigHig[key] = np.asarray(exp_2sigHig[key], dtype='f')
            for im, m in enumerate(masses_[iMVAmass]):
                graph_obs[key].SetPoint(im, np.asarray([m], dtype='f')[0], obs[key][im])
                graph_exp[key].SetPoint(im, np.asarray([m], dtype='f')[0], exp[key][im])
#                print key, im, m, exp[key][im]
                graph_exp1sig[key].SetPoint(im, np.asarray([m], dtype='f')[0], exp_1sigLow[key][im])
                graph_exp1sig[key].SetPoint(2* len(masses_[iMVAmass]) -1 -im, np.asarray([m], dtype='f')[0], exp_1sigHig[key][im])
                graph_exp2sig[key].SetPoint(im, np.asarray([m], dtype='f')[0], exp_2sigLow[key][im])
                graph_exp2sig[key].SetPoint(2* len(masses_[iMVAmass]) -1 -im, np.asarray([m], dtype='f')[0], exp_2sigHig[key][im])
        
                # Latex dump
                if latex:
                    expected = exp[key][im]
                    with open('tables/%s.tex' % (outfilename), 'a') as f:
                        f.write( (R"%d & $%.2e}$ & $%.2e}$ & $-%.2e}$ / $+%.2e}$ & $-%.2e}$ / $+%.2e}$ \\" % (m, obs[key][im], expected, expected - exp_1sigLow[key][im], exp_1sigHig[key][im] - expected, expected - exp_2sigLow[key][im], exp_2sigHig[key][im] - expected) + '\n').replace('e+0', R'\times 10^{') )
                if dumpForCommonPlot:
                    with open(dumpForCommonPlotOutfile, 'a') as dumpfile:
                        dumpfile.write('%i\t%.1f\t%.1f\n' % (m, obs[key][im], exp[key][im]))

            graph_obs_fullsim[key] = TGraph()
            ipoint = 0
            for im, m in enumerate(masses_[iMVAmass]):
                if m not in masses[iMVAmass]:
                    continue
                graph_obs_fullsim[key].SetPoint(ipoint, np.asarray([m], dtype='f')[0], obs[key][im])
                ipoint += 1
                
            
            graph_obs[key].SetLineWidth(2)
            graph_obs[key].SetLineStyle(1)
            graph_obs[key].SetMarkerStyle(1)
            graph_obs_fullsim[key].SetLineWidth(2)
            graph_obs_fullsim[key].SetLineStyle(1)
            graph_obs_fullsim[key].SetMarkerSize(0.8)
            graph_obs_fullsim[key].SetMarkerStyle(20)
            graph_exp[key].SetLineWidth(2)
            graph_exp[key].SetLineStyle(3)
            if "cut-n-count" in key:
                graph_obs[key].SetLineColor(ROOT.kBlack)
                graph_obs_fullsim[key].SetLineColor(ROOT.kBlack)
                graph_exp[key].SetLineColor(ROOT.kBlack)
                graph_exp1sig[key].SetFillColor(color.GetColor('#77FF33'))
                graph_exp1sig[key].SetLineColor(color.GetColor('#77FF33'))
                graph_exp2sig[key].SetFillColor(ROOT.kYellow-4)
                graph_exp2sig[key].SetLineColor(ROOT.kYellow-4)
            elif key == "shape-hist":
                graph_exp[key].SetLineColor(ROOT.kMagenta+2)
            graph_exp[key].GetXaxis().SetLimits(xmin,xmax)
            graph_exp[key].SetMaximum(ymax)
            graph_exp[key].SetMinimum(ymin)
            graph_exp[key].GetXaxis().SetTitle("m_{X}^{spin-%i}(GeV)" % spin)
            graph_exp[key].GetYaxis().SetTitleSize(0.032)
            if not scaledToHH:
                yaxisextra = " #times BR(HH #rightarrow b#bar{b}l#nul#nu)"
            graph_exp[key].GetYaxis().SetTitle("95% CL limit on #sigma(pp#rightarrow X^{spin-" + str(spin) + "} #rightarrow HH)" + yaxisextra + " (fb)")

    graph_th = {}
    xsec_x_br = {}
    xsec_x_br_params = {}
    for theoryFile in theoryCurves[ispin]:
        graph_th[theoryFile] = TGraph()
        print '# Crossings with theory= %s' % theoryFile
        with open(theoryFile) as f:
            alldata = json.load(f)
            xsec_x_br[theoryFile] = alldata['xsec_x_br']
            xsec_x_br_params[theoryFile] = alldata['parameters']
        l_m = [x for x in xsec_x_br[theoryFile]]
        l_m = map(int, l_m)
        l_m.sort()
        exp_minus_th = 0
        obs_minus_th = 0
        for im, m in enumerate(l_m):
            y = xsec_x_br[theoryFile][unicode(m)]
            if not scaledToHH and not xsec_x_br_params[theoryFile]['scaled_to_final_state']:
                y *= BRbbVV
            elif scaledToHH and xsec_x_br_params[theoryFile]['scaled_to_final_state']:
                y /= BRbbVV
            if fbUnits and 'pb' in str(xsec_x_br_params[theoryFile]['unit']):
                y *= 1.0e3
            elif not fbUnits and 'fb' in str(xsec_x_br_params[theoryFile]['unit']):
                y /= 1.0e3
            graph_th[theoryFile].SetPoint(im, m, y)
            exp_ = 0
            obs_ = 0
            if im == 0:
                for iMVAmass, MVAmass in enumerate(MVAmasses):
                    for ikey, key in enumerate(strategies):
                        if str(MVAmass) not in key:
                            continue
                        if m in masses_[iMVAmass]:
                            imass = masses_[iMVAmass].index(m)
                            obs_ = obs[key][imass]
                            exp_ = exp[key][imass]
                exp_minus_th = exp_ - y
                obs_minus_th = obs_ - y
            else:
                for iMVAmass, MVAmass in enumerate(MVAmasses):
                    for ikey, key in enumerate(strategies):
                        if str(MVAmass) not in key:
                            continue
                        if m in masses_[iMVAmass]:
                            imass = masses_[iMVAmass].index(m)
                            obs_ = obs[key][imass]
                            exp_ = exp[key][imass]
                # if there is a crossing, just dump the value on screen
                if 260 < m <= 900 and (exp_ - y) * exp_minus_th < 0:
                    print '# \tAt m= %i, crossing with exp= %f (obs= %f)' % (m, exp_, obs_)
                if 260 < m <= 900 and (obs_ - y) * obs_minus_th < 0:
                    print '# \tAt m= %i, crossing with obs= %f (exp= %f)' % (m, obs_, exp_)
                exp_minus_th = exp_ - y
                obs_minus_th = obs_ - y
                
            

    # First uncertainty bands
    for ikey, key in enumerate(strategies):
        if ikey == 0:
            graph_exp[key].Draw("al")
        else:
            graph_exp[key].Draw("l")
        graph_exp2sig[key].Draw("f")
        graph_exp1sig[key].Draw("f")
    # then theory curves
    for ith, theoryFile in enumerate(theoryCurves[ispin]):
        graph_th[theoryFile].SetLineWidth(2)
        # colors from http://www.colourlovers.com/palette/694737/Thought_Provoking 
        if ith == 0:
#            graph_th[theoryFile].SetLineColor(color.GetColor('#ECD078'))
            graph_th[theoryFile].SetLineColor(color.GetColor('#FF5533'))
        elif ith == 1:
            graph_th[theoryFile].SetLineColor(color.GetColor('#3388FF'))
        elif ith == 2:
            graph_th[theoryFile].SetLineColor(color.GetColor('#FF3399'))
        elif ith == 3:
            graph_th[theoryFile].SetLineColor(color.GetColor('#53777A'))
        graph_th[theoryFile].SetLineStyle(5)
        graph_th[theoryFile].Draw('l')
    # then median expected and observed
    for ikey, key in enumerate(strategies):
        graph_exp[key].Draw("l")
        if unblind:
            graph_obs[key].Draw("l")
            graph_obs_fullsim[key].Draw('p')
    for xline in transitionLines:
        line = TLine(xline, ymin, xline, ymax)
        line.SetLineStyle(7)
        line.Draw('same')

    legend = TLegend(0.47,0.55,0.90,0.92)
    legend.SetTextFont(42)
    if scaledToHH:
        legend.SetHeader('X^{spin-' + str(spin) + '} #rightarrow HH #rightarrow b#bar{b}l#nul#nu analysis')
    #legend.SetTextSize(0.028)
    legend.SetFillStyle(0)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetLineColor(ROOT.kWhite)
    legend.SetShadowColor(ROOT.kWhite)
    legend.AddEntry(graph_obs_fullsim[strategies[0]], "Observed 95% upper limit", "lp")
    legend.AddEntry(graph_exp[strategies[0]], "Expected 95% upper limit", "l")
    legend.AddEntry(graph_exp1sig[strategies[0]], "Expected limit #pm 1 std. deviation", "f")
    legend.AddEntry(graph_exp2sig[strategies[0]], "Expected limit #pm 2 std. deviations", "f")
    first2HDM = False
    for theoryFile in theoryCurves[ispin]:
        th_tlatex = xsec_x_br_params[theoryFile]['TLatex']
        if len(th_tlatex) == 1:
            legend.AddEntry(graph_th[theoryFile], '', '')
            legend.AddEntry(graph_th[theoryFile], th_tlatex[0], 'l')
        elif len(th_tlatex) == 2 and '2HDM' in th_tlatex[0] and not first2HDM:
            legend.AddEntry(graph_th[theoryFile], '', '')
            legend.AddEntry(graph_th[theoryFile], '#splitline{%s}{%s}' % (th_tlatex[0].replace(' Type I', ''), th_tlatex[1]), '')
            first2HDM = True
        if len(th_tlatex) == 2 and '2HDM' in th_tlatex[0]:
            legend.AddEntry(graph_th[theoryFile], th_tlatex[0].split(',')[0].split('2HDM ')[1], 'l')
    legend.Draw("same")
    
    # more cosmetics
    intL = 2245.792 # in pb-1
    latexLabel = TLatex()
    latexLabel.SetTextSize(0.75 * c1.GetTopMargin())
    latexLabel.SetNDC()
    latexLabel.SetTextFont(42) # helvetica
    latexLabel.DrawLatex(0.69, 0.96, "2.30 fb^{-1} (13 TeV)")
    latexLabel.SetTextFont(61) # helvetica bold face
    latexLabel.DrawLatex(0.13, 0.96, "CMS")
    latexLabel.SetTextFont(52) # helvetica italics
    latexLabel.DrawLatex(0.22, 0.96, "Preliminary")
    if scaledToHH:
        latexLabel.DrawLatex(0.16, 0.88, "Assumes SM BR")
    
    c1.Print("figures/%s.png" % (outfilename))
    c1.Print("figures/%s.pdf" % (outfilename))
    
    with open('tables/%s.tex' % (outfilename), 'a') as f:
        f.write(R'\bottomrule' + '\n')
        f.write(R'\end{tabular}' + '\n')
    


