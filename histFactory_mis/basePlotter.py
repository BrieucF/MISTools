import copy, sys
    #from HHAnalysis import HH
    #from ScaleFactors import *


class BasePlotter:
#    def __init__(self):

    def generatePlots(self, categories = ["All"], stage = "no_cut", requested_plots = ["bdtoutput"], weightExpression = "", extraCut = "", systematic = "nominal", extraString = ""):

        # MVA evaluation : ugly but necessary part
        baseStringForMVA_part1 = 'evaluateMVA("/home/fynu/bfrancois/scratch/llbb_release/modifieForJetEtaTooHigh/CMSSW_5_3_20/src/UserCode/ttbar_effth_delphes/analyzer/weights/NAME.weights.xml", '
        MVAsuffixs = ["_vs_All_relativeWeight_withZh_kBDT", "_vs_All_weight1_withZh_kBDT"]
        MVAsuffixs_noZh = ["_vs_All_relativeWeight_withoutZh_kBDT", "_vs_All_weight1_withoutZh_kBDT"]
        baseStringForMVA_part2 = '{{"MinusLogW_gg_Zbb", MinusLogW_gg_Zbb}, {"MinusLogW_TT", MinusLogW_TT}, {"MinusLogW_ZZ_cor0", MinusLogW_ZZ_cor0}, {"MinusLogW_ZH_cor0", MinusLogW_ZH_cor0}})'
        baseStringForMVA_part2_noZh = '{{"MinusLogW_gg_Zbb", MinusLogW_gg_Zbb}, {"MinusLogW_TT", MinusLogW_TT}, {"MinusLogW_ZZ_cor0", MinusLogW_ZZ_cor0}})'
        stringForMVA = baseStringForMVA_part1 + baseStringForMVA_part2
        BDToutputs = {}
        bdtNames = []
        BDToutputsVariable = {}
        processes = ["DY","TT","ZZ","ZH"]
        processes_noZh = ["DY","TT","ZZ"]
        for suffix in MVAsuffixs_noZh :
            for proc in processes_noZh :
                bdtName = "BDT_" + proc + suffix
                bdtNames.append(bdtName)
                BDToutputsVariable[bdtName] = baseStringForMVA_part1.replace("NAME", bdtName) + baseStringForMVA_part2_noZh
        for suffix in MVAsuffixs :
            for proc in processes :
                bdtName = "BDT_" + proc + suffix 
                bdtNames.append(bdtName)
                BDToutputsVariable[bdtName] = baseStringForMVA_part1.replace("NAME", bdtName) + baseStringForMVA_part2

        summed_withZhBDTPlotNames = []
        for suffix in MVAsuffixs :
            iterator = 0 
            plotName = "summed"+suffix.replace("vs_All_","")
            summed_withZhBDTPlotNames.append(plotName)
            BDToutputsVariable[plotName] = ""
            for proc in processes :
                iterator += 1
                bdtName = "BDT_" + proc + suffix
                if iterator == len(processes) :
                    BDToutputsVariable[plotName] += baseStringForMVA_part1.replace("NAME", bdtName) + baseStringForMVA_part2
                else : 
                    BDToutputsVariable[plotName] += baseStringForMVA_part1.replace("NAME", bdtName) + baseStringForMVA_part2 + " + "
                    
        multiplied_withZhBDTPlotNames = []
        for suffix in MVAsuffixs :
            iterator = 0 
            plotName = "multiplied"+suffix.replace("vs_All_","")
            multiplied_withZhBDTPlotNames.append(plotName)
            BDToutputsVariable[plotName] = ""
            for proc in processes :
                iterator += 1
                bdtName = "BDT_" + proc + suffix
                if iterator == len(processes) :
                    BDToutputsVariable[plotName] += baseStringForMVA_part1.replace("NAME", bdtName) + baseStringForMVA_part2
                else : 
                    BDToutputsVariable[plotName] += baseStringForMVA_part1.replace("NAME", bdtName) + baseStringForMVA_part2 + " * "

        summed_withoutZhBDTPlotNames = []
        for suffix in MVAsuffixs_noZh :
            iterator = 0 
            plotName = "summed"+suffix.replace("vs_All_","")
            summed_withoutZhBDTPlotNames.append(plotName)
            BDToutputsVariable[plotName] = ""
            for proc in processes_noZh :
                iterator += 1
                bdtName = "BDT_" + proc + suffix
                if iterator == len(processes_noZh) :
                    BDToutputsVariable[plotName] += baseStringForMVA_part1.replace("NAME", bdtName) + baseStringForMVA_part2_noZh
                else : 
                    BDToutputsVariable[plotName] += baseStringForMVA_part1.replace("NAME", bdtName) + baseStringForMVA_part2_noZh + " + "
                    
        multiplied_withoutZhBDTPlotNames = []
        for suffix in MVAsuffixs_noZh :
            iterator = 0 
            plotName = "multiplied"+suffix.replace("vs_All_","")
            multiplied_withoutZhBDTPlotNames.append(plotName)
            BDToutputsVariable[plotName] = ""
            for proc in processes_noZh :
                iterator += 1
                bdtName = "BDT_" + proc + suffix
                if iterator == len(processes_noZh) :
                    BDToutputsVariable[plotName] += baseStringForMVA_part1.replace("NAME", bdtName) + baseStringForMVA_part2_noZh
                else : 
                    BDToutputsVariable[plotName] += baseStringForMVA_part1.replace("NAME", bdtName) + baseStringForMVA_part2_noZh + " * "


                

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
        self.dict_cat_cut["All"] = cut_for_All_channel

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
                        'binning' : '(60, -1, 1)'
                })
            for plotName in summed_withZhBDTPlotNames :
                self.bdtoutput_plot.append({
                        'name' : 'MVA_%s_%s_%s%s'%(plotName, self.llFlav, self.extraString, self.systematicString),
                        'variable' : BDToutputsVariable[plotName],
                        'plot_cut' : self.totalCut,
                        'binning' : '(60, -0.5, 1.5)'
                })
            for plotName in multiplied_withZhBDTPlotNames :
                self.bdtoutput_plot.append({
                        'name' : 'MVA_%s_%s_%s%s'%(plotName, self.llFlav, self.extraString, self.systematicString),
                        'variable' : BDToutputsVariable[plotName],
                        'plot_cut' : self.totalCut,
                        'binning' : '(60, -0.1, 0.1)'
                })
            for plotName in summed_withoutZhBDTPlotNames :
                self.bdtoutput_plot.append({
                        'name' : 'MVA_%s_%s_%s%s'%(plotName, self.llFlav, self.extraString, self.systematicString),
                        'variable' : BDToutputsVariable[plotName],
                        'plot_cut' : self.totalCut,
                        'binning' : '(60, -1, 1)'
                })
            for plotName in multiplied_withoutZhBDTPlotNames :
                self.bdtoutput_plot.append({
                        'name' : 'MVA_%s_%s_%s%s'%(plotName, self.llFlav, self.extraString, self.systematicString),
                        'variable' : BDToutputsVariable[plotName],
                        'plot_cut' : self.totalCut,
                        'binning' : '(60, -0.2, 0.2)'
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

