import copy, sys
#from HHAnalysis import HH
#from ScaleFactors import *


class BasePlotter:
#    def __init__(self):

    def generatePlots(self, categories = ["All"], stage = "no_cut", requested_plots = ["bdtoutput"], weightExpression = "(runOnMC)? (btaggingReweightingMM*leptonsReweightingweight*lumiReweightingLumiWeight*mcReweightingweight) : (1.)", extraCut = "", systematic = "nominal", extraString = ""):

        # MVA evaluation : ugly but necessary part
        baseStringForMVA_part1 = 'evaluateMVA("/home/fynu/bfrancois/scratch/llbb_release/modifieForJetEtaTooHigh/CMSSW_5_3_20/src/UserCode/ttbar_effth_delphes/analyzer/weights/NAME.weights.xml", '
        MVAsuffixs = ["_vs_All_relativeWeight_withZh_kBDT"]
        baseStringForMVA_part2 = '{{"MinusLogW_gg_Zbb", MinusLogW_gg_Zbb}, {"MinusLogW_TT", MinusLogW_TT}, {"MinusLogW_ZZ_cor0", MinusLogW_ZZ_cor0}, {"MinusLogW_ZH_cor0", MinusLogW_ZH_cor0}})'
        stringForMVA = baseStringForMVA_part1 + baseStringForMVA_part2
        BDToutputs = {}
        bdtNames = []
        BDToutputsVariable = {}
        processes = ["DY","TT","ZZ","ZH"]
        for proc in processes :
            for suffix in MVAsuffixs :
                bdtName = "BDT_" + proc + suffix 
                bdtNames.append(bdtName)
                BDToutputsVariable[bdtName] = baseStringForMVA_part1.replace("NAME", bdtName) + baseStringForMVA_part2

        # Possible stages
        dict_stage_cut = {
               "no_cut" : "", 
               }
        # Categories
        self.dict_cat_cut =  {
            "ElEl" : "(IncLepslep1isMuon && IncLepslep2isMuon)",
            "MuMu" : "(IncLepslep1isElectron && IncLepslep2isElectron)",
            "MuEl" : "((IncLepslep1isElectron && IncLepslep2isMuon) || (IncLepslep1isMuon && IncLepslep2isElectron))"
            }   
        cut_for_All_channel = "(" + self.dict_cat_cut["ElEl"] + "||" + self.dict_cat_cut["MuMu"] + "||" +self.dict_cat_cut["MuEl"] + ")"
        cut_for_SF_channel = "(" + self.dict_cat_cut["ElEl"] + "||" + self.dict_cat_cut["MuMu"] + ")"
        self.dict_cat_cut["SF"] = cut_for_SF_channel
        self.dict_cat_cut["All"] = "rc_stage_0_idx" #cut_for_All_channel

        ###########
        # Weights #
        ###########

        # Append the proper extension to the name plot if needed (scale name are at the end of the code)
        self.systematicString = ""
        if not systematic == "nominal" and not systematic == "scale" :
            self.systematicString = "__" + systematic

        #########
        # PLOTS #
        #########
        self.bdtoutput_plot = []
        for cat in categories :

            catCut = self.dict_cat_cut[cat]
            self.totalCut = self.joinCuts(catCut, extraCut, dict_stage_cut[stage])
            self.llFlav = cat
            self.extraString = stage + extraString

            # BDT output plots
            for bdtName in bdtNames :
                self.bdtoutput_plot.append({
                        'name' : 'MVA_%s_%s_%s%s'%(bdtName, self.llFlav, self.extraString, self.systematicString),
                        'variable' : BDToutputsVariable[bdtName],
                        'plot_cut' : self.totalCut,
                        'binning' : '(60, -0.7, 0.7)'
                })

        plotsToReturn = []
        for plotFamily in requested_plots :
            for plot in getattr(self, plotFamily+"_plot"):
                plot["weight"] = weightExpression 
                plotsToReturn.append(plot)

        return plotsToReturn


    def joinCuts(self, *cuts):
        if len(cuts) == 0: 
            return ""
        elif len(cuts) == 1: 
            return cuts[0]
        else:
            totalCut = "("
            for cut in cuts:
                cut = cut.strip().strip("&")
                if cut == "":
                    continue
                totalCut += "(" + cut + ")&&" 
            totalCut = totalCut.strip("&") + ")"
            return totalCut

