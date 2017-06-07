import sys
import ROOT
import numpy as np
import pprint

rooFileName = sys.argv[1]
print rooFileName
rootFile = ROOT.TFile(rooFileName)
tree_fit_b = rootFile.Get("tree_fit_b")
proc_yield = {"DY" : 2188.74, "TT" : 11123.8, "tWp" : 194.418, "tWm" : 190.568, "ZZ" : 49.4189, "ZH" : 26.8249}
proc_yields_toys = {} # Mean, std deviation, ratio, relative std deviation
for proc in proc_yield.keys() :
    list_yields = []
    count = 0.0
    liefName = "n_exp_binyields_proc_" + proc
    if not liefName in tree_fit_b.GetListOfBranches() : 
        print "No lief name : ", liefName
        liefName = "n_exp_final_binyields_proc_" + proc
    for entry in tree_fit_b :
        if entry.fit_status >= 0 :
            count += 1
            list_yields.append(getattr(entry, liefName))

    proc_yields_toys[proc] = [np.mean(list_yields), np.std(list_yields), np.mean(list_yields)/proc_yield[proc],  np.std(list_yields)/np.mean(list_yields)]

pp = pprint.PrettyPrinter(indent=4)
pp.pprint(proc_yields_toys)
print count, " fit succeded (%s %%)"%(count*100/float(tree_fit_b.GetEntries()))
print "Sum of absolute errors : ", sum([i for i in [proc_yields_toys[j][1] for j in proc_yields_toys.keys()]])
print "Sqrt(Sum of (absolute errors)^2) : ", np.sqrt(sum([i*i for i in [proc_yields_toys[j][1] for j in proc_yields_toys.keys()]]))
print "Sum of relative errors : ", sum([i for i in [proc_yields_toys[j][3] for j in proc_yields_toys.keys()]])
print "Sqrt(Sum of (relative errors)^2) : ", np.sqrt(sum([i*i for i in [proc_yields_toys[j][3] for j in proc_yields_toys.keys()]]))


