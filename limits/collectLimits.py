#! /bin/env python

import os, sys, argparse, json
import stat
import subprocess
import glob
import re

import numpy as np

import ROOT
ROOT.gROOT.SetBatch()
ROOT.PyConfig.IgnoreCommandLineOptions = True

# Open /dev/null to redirect stdout
DEVNULL = open(os.devnull, 'wb')

def getLimitsFromFile(input_file):
    """
    Extract observed, expected, and 1/2 sigma limits from an HybridNew grid
    """

    f = ROOT.TFile.Open(input_file)
    if not f or f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered):
        return None

    data = {}
    # Index 0 is DOWN error, index 1 is UP error
    one_sigma = data.setdefault('one_sigma', [0, 0])
    two_sigma = data.setdefault('two_sigma', [0, 0])

    limit = f.Get('limit')

    limit.GetEntry(2)
    data['expected'] = round(limit.limit, 2)

    limit.GetEntry(5)
    data['observed'] = round(limit.limit, 2)

    limit.GetEntry(4)
    two_sigma[1] = round(limit.limit - data['expected'], 2)

    limit.GetEntry(0)
    two_sigma[0] = round(data['expected'] - limit.limit, 2)

    limit.GetEntry(3)
    one_sigma[1] = round(limit.limit - data['expected'], 2)

    limit.GetEntry(1)
    one_sigma[0] = round(data['expected'] - limit.limit, 2)

    return data

def transform_param(p):
    return float(p.replace("p", ".").replace("m", "-"))

# Options
parser = argparse.ArgumentParser(description='Collection non-resonant limits')
parser.add_argument('-d', '--directory', action='store', type=str, default="./", help='Directory where limits are')
parser.add_argument('-o',  '--output', action='store', type=str, default='limits.json', help='Name of output file (JSON format)')

options = parser.parse_args()

print("Extracting limits...")

# Extract limits for all points, and save into a JSON file
limits = []
limits_dict_mis = {}

m = re.compile('higgsCombineMIS_Signal_.+.Asymptotic')
for directory in os.listdir(options.directory):
    path = os.path.join(options.directory, directory)
    if os.path.isdir(path):
        for file in os.listdir(path):
            if m.search(file):
                print "Treating ",  file
                point_limits = getLimitsFromFile(os.path.join(path,file))


                if point_limits['expected'] == 0:
                    print("Warning: expected is 0, skipping point")
                    continue
                limits_dict_mis[directory] = point_limits
                limits.append({
                    'signal': directory,
                    'limits': point_limits
                    })
limits_dict_dedicated = {}
limits_dict_dedicated["TTDM_1_10"] = {"expected": 25467}
limits_dict_dedicated["TTDM_1_100"] = {"expected": 9412}
limits_dict_dedicated["TTDM_1_500"] = {"expected": 2051}
limits_dict_dedicated["TTDM_1_10_ps"] = {"expected": 10141}
limits_dict_dedicated["TTDM_1_100_ps"] = {"expected": 4200}
limits_dict_dedicated["TTDM_1_500_ps"] = {"expected": 2271}
limits_dict_dedicated["X_400_hh"] = {"expected": 259}
limits_dict_dedicated["X_650_hh"] = {"expected": 90.75}
limits_dict_dedicated["X_900_hh"] = {"expected": 75.75}
limits_dict_dedicated["SMhh"] = {"expected": 92.8}
limits_dict_dedicated["H_500_ZA_300"] = {"expected": 237}
limits_dict_dedicated["H_800_ZA_700"] = {"expected": 116.7}
limits_dict_dedicated["Stop_500_325"] = {"expected": 3453.1}
limits_dict_dedicated["Stop_850_100"] = {"expected": 2335.9}

output_file = options.output
with open(output_file, 'w') as f:
    json.dump(limits, f)
print("Limits saved as %s" % output_file)
#Print it as a latex table!
print r"""
\begin{table}
\label{tab:limits_mis}
\begin{center}
\begin{tabular}{@{}lcc@{}} \toprule
Signal model (Masses in GeV)&CLs (dedicated search)& CLs (MIS)\\ \midrule
$t\bar t$ + scalar DM (m$_\chi$=1, m$_\phi$=10)  & $25467$ & $%s$  \\ 
$t\bar t$ + scalar DM (m$_\chi$=1, m$_\phi$=100)  & $9412$ & $%s$  \\ 
$t\bar t$ + scalar DM (m$_\chi$=1, m$_\phi$=500)  & $2051$ & $%s$  \\ 
%%$t\bar t$ + pseudo scalar DM (m$_\chi$=1, m$_\phi$=10)  & $10141$ & $%s$  \\ 
%%$t\bar t$ + pseudo scalar DM (m$_\chi$=1, m$_\phi$=100)  & $4200$ & $%s$  \\ 
%%$t\bar t$ + pseudo scalar DM (m$_\chi$=1, m$_\phi$=500)  & $2271$ & $%s$  \\ 
\midrule
$X \to hh$  m$^{ spin-0}_{X} = 400$ & $259$ & $%s$  \\ 
$X \to hh$ m$^{ spin-0}_X = 650$ & $90.75$ & $%s$  \\ 
$X \to hh$ m$^{ spin-0}_X = 900$ & $75.75$ & $%s$  \\ 
\midrule
SM $hh$ & $92.8$ & $%s$ \\
\midrule
$H(500) \to Z ~ A(300) $ & $237$ & $%s$ \\
$H(800) \to Z ~ A(700) $ & $116.7$ & $%s$ \\
\midrule
$\tilde{t}\tilde{t}$ $(m_{\tilde t}=500, m_{\chi} = 325)$ & 3453.1 & $%s$ \\
$\tilde{t}\tilde{t}$ $(m_{\tilde t}=850, m_{\chi} = 100)$ &  2335.9 & $%s$ \\
\bottomrule
\end{tabular}
\caption{95\%% expected CLs limits in $fb$ obtained by the dedicated searches (middle raw) and by the MIS method using the yields in each of the final boxes. All the systematic uncertainties described in Sec.~\ref{sec:syst} are taken into account except for the additional global scaling uncertainties. The cross section uncertainty on \ttbar and DY are replaced by the uncertainty on the background normalization fit. Each line in the table corresponds to the limit on a different BSM scenarios/parameter.}
\end{center}
\end{table}
"""%(   limits_dict_mis["TTDM_1_10"]["expected"], 
        limits_dict_mis["TTDM_1_100"]["expected"], 
        limits_dict_mis["TTDM_1_500"]["expected"], 
        limits_dict_mis["TTDM_1_10_ps"]["expected"], 
        limits_dict_mis["TTDM_1_100_ps"]["expected"], 
        limits_dict_mis["TTDM_1_500_ps"]["expected"], 
        limits_dict_mis["X_400_hh"]["expected"], 
        limits_dict_mis["X_650_hh"]["expected"], 
        limits_dict_mis["X_900_hh"]["expected"], 
        limits_dict_mis["SMhh"]["expected"], 
        limits_dict_mis["H_500_ZA_300"]["expected"], 
        limits_dict_mis["H_800_ZA_700"]["expected"], 
        limits_dict_mis["Stop_500_325"]["expected"], 
        limits_dict_mis["Stop_850_100"]["expected"],
    )
ratio = []
ratio_dict = {}
for signal in  limits_dict_dedicated.keys():
    #if not "Stop" in signal:
        ratio_dict[signal] = limits_dict_mis[signal]["expected"]/limits_dict_dedicated[signal]["expected"]
        ratio.append(limits_dict_mis[signal]["expected"]/limits_dict_dedicated[signal]["expected"])
print ratio
print ratio_dict
print min(ratio)
print max(ratio)

#for i in options.inputFolder:
#
#    print i
#    res = m.search(i)
#    if not res:
#        print("Warning: ignoring input {}. Unexpected filename format".format(i))
#        continue
#
#    kl = res.group(1)
#    kt = res.group(2)
#
#    kl = transform_param(kl)
#    kt = transform_param(kt)
#
#    point_limits = getLimitsFromFile(i)
#
#    if point_limits['expected'] == 0:
#        print("Warning: expected is 0, skipping point")
#        continue
#
#    limits.append({
#        'parameters': (kl, kt),
#        'limits': point_limits
#        })
#
#
#output_file = options.output
#with open(output_file, 'w') as f:
#    json.dump(limits, f)
#
#print("Limits saved as %s" % output_file)
#
