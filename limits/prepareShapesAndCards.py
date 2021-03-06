#! /bin/env python

# Python imports
import os, sys, stat, argparse, getpass, json
from datetime import datetime
from math import sqrt

# to prevent pyroot to hijack argparse we need to go around
tmpargv = sys.argv[:] 
sys.argv = []

# ROOT imports
import ROOT
ROOT.gROOT.SetBatch()
ROOT.PyConfig.IgnoreCommandLineOptions = True
sys.argv = tmpargv

def zeroNegativeBins(h):
    for i in range(1, h.GetNbinsX() + 1):
        if h.GetBinContent(i) < 0.:
            h.SetBinContent(i, 0.)

def get_options():
    """
    Parse and return arguments provided by the user
    """
    parser = argparse.ArgumentParser(description='Create shape datacards ready for combine')

    #parser.add_argument('-p', '--path', action='store', dest='root_path', type=str, default='/home/fynu/bfrancois/scratch/framework/MIS_prod_data/CMSSW_7_6_5/src/cp3_llbb/HHTools/histFactory_hh/weightPlot_v4/condor/output/', help='histFactory output path')
    parser.add_argument('-p', '--path', action='store', dest='root_path', type=str, default='/home/fynu/bfrancois/scratch/framework/MIS_prod_data/CMSSW_7_6_5/src/cp3_llbb/ttbar_effth_delphes/analyzer/results/one_vs_all_puregini_10pc_effratio_2_nonInfPunderatedWeight_DYdfWeight/suitable_th1_ordered', help='histFactory output path')

    LUMI = 2300
    parser.add_argument('-l', '--luminosity', action='store', type=float, dest='luminosity', default=LUMI, help='integrated luminosity (default is %f /pb)' % LUMI)

    # hypothesis options
    #parser.add_argument('-c', '--cluster', nargs='+', metavar='CLUSTER', type=str, dest='clusters', default=['SM'], help='clusters hypothesis. Use `all` as an alias for all the clusters')

    # MVA options
    #parser.add_argument('-m', '--trained-cluster', nargs='+', type=str, dest='MVAclusters',
    #    default=['SM'], help='MVA training clusters. Use `all` as an alias for all MVA clusters, or `analysis` to use analysis conditions')

    parser.add_argument('-o', '--output', action='store', dest='output', type=str, default='Cards_yields', help='Output directory')
    #parser.add_argument('-n', '--name', action='store', dest='name', type=str, default='dummy', help='Name of the box splitting, only relevant for box mode')
    #parser.add_argument('-f', '--file', action='store', dest='file', type=str, default='', help='Name of the input rootfile, relevant only for box mode')
    parser.add_argument('-b', '--box', action='store_true', dest='box', help='Extract limits based on box yields')
    parser.add_argument('--reweight', action='store_true', dest='reweight', help='Apply a preliminary reweighting')

    parser.add_argument('--fake-data', action='store_true', dest='fake_data',
        help='Use fake data instead of real data')
    parser.add_argument('--SF', action='store_true', dest='SF',
        help='Produce cards for scale factors extraction (add line with rateParam)')

    parser.add_argument('--nosys', action='store', dest='nosys', default=False,
            help='Consider or not systematic uncertainties (NB : bbb uncertainty is with another flag)')
    parser.add_argument('--bbb', action='store', dest='bbb', default=False,
        help='Consider or not bin by bin MC stat systematic uncertainties')
    parser.add_argument('--mt', action='store', dest='mt', type=float, default=0.5,
        help='Consider or not bin by bin MC stat systematic uncertainties')
    parser.add_argument('--flavSplit', action='store', dest='flavSplit', type=bool, default=False, 
            help='Split flavor or not')

    options = parser.parse_args()
    return options

def main():
    """Main function"""
    global options

    # get the options
    options = get_options()

    if options.box and options.SF :
        backgrounds = ['ZZ', 'TT', 'DY', 'ZH', 'tWp', 'tWm']
        #backgrounds = ['ZZ', 'TT', 'DY', 'DY_10_50', 'ZH', 'tWp', 'tWm', 'TTV', 'ST', 'WW', 'ttH']
        signals = [
                    #'TTDM_1_10',
                    #'TTDM_1_100',
                    #'TTDM_1_500',
                    #'TTDM_1_10_ps',
                    #'TTDM_1_100_ps',
                    #'TTDM_1_500_ps',
                    #'H_200_ZA_50',
                    #'H_3000_ZA_2000',
                    #'H_500_ZA_300',
                    #'H_800_ZA_700',
                    #'Stop_500_325',
                    #'Stop_850_100',
                    #'X_400_hh',
                    #'X_650_hh',
                    #'X_900_hh',
                    'SMhh',
                    ]
    else :
        #backgrounds = ['ZZ', 'TT', 'DYbb', 'DYbx', 'DYxx', 'DY_10_50', 'ZH', 'tWp', 'tWm', 'TTV', 'ST', 'WW', 'ttH']
        backgrounds = ['ZZ', 'TT', 'DY', 'DY_10_50', 'ZH', 'tWp', 'tWm', 'TTV', 'ST', 'WW', 'ttH']
        #backgrounds = ['ZZ', 'TT', 'DY', 'ZH', 'tWp', 'tWm']
        signals = [
                    'TTDM_1_10',
                    'TTDM_1_100',
                    'TTDM_1_500',
                    'TTDM_1_10_ps',
                    'TTDM_1_100_ps',
                    'TTDM_1_500_ps',
                    #'H_200_ZA_50',
                    #'H_3000_ZA_2000',
                    'H_500_ZA_300',
                    'H_800_ZA_700',
                    'Stop_500_325',
                    'Stop_850_100',
                    'X_400_hh',
                    'X_650_hh',
                    'X_900_hh',
                    'SMhh',
                    ]

    if options.box :
        print "Box mode!"
        if options.flavSplit :
            discriminants = { 
                     #'yields': [(1, 'yields')],
                     'yields_flavSplit': [(1, 'yields_sf'), (2, 'yields_df')],
                     #'ordered_yields': [(1, 'ordered_yields')]
                     }
        else :
            discriminants = {
                'yields': [(1, 'yields')],
                }
        for discriminant in discriminants.keys() :
            prepareShapes(backgrounds, signals, discriminants[discriminant], discriminant)
        #prepareShapesFromBoxes(backgrounds, signals, options.name)

    else :
        print "Variable mode!"
        discriminants = { 
                 #"ATan_DY_TT" : [(1, 'ATan_DY_TT')],
                 #"mll" : [(1, 'mll')],
                 #"mll_ee" : [(1, 'mll_ee')],
                 #"mll_mumu" : [(1, 'mll_mumu')],
                 "mll_emu" : [(1, 'mll_emu')],
                 #"ATan_DY_TT_plus_CSVProd" : [(1, 'ATan_DY_TT'), (2, 'csvProd')],
                 #"CSVProd" : [(1, 'csvProd')],
                 #"ATan_twp_twm" : [(1, 'ATan_twp_twm')],
                 #"tt_dy_weights": [(1, 'TT_weights'), (2, 'DY_weights')],
                 #"dy_weights": [(2, 'DY_weights')],
                 #"bare_weights": [(1, 'TT_weights'), (2, 'DY_weights'), (3, 'ZH_weights'), (4, 'tWp_weights'), (5, 'tWm_weights'), (6, 'ZZ_weights')],

                 #"ThreeBDTTogether_noZh" : [(1, 'TT_vs_All_noZh'), (2, 'DY_vs_All_noZh'), (3, 'ZZ_vs_All_noZh')],
                 #"ThreeBDTTogether_Zh" : [(1, 'TT_vs_All_Zh'), (2, 'DY_vs_All_Zh'), (3, 'ZZ_vs_All_Zh'), (4, 'ZH_vs_All_Zh')],
                 #"BDTsummed_noZh" : [(1, 'BDTsummed_noZh')],
                 #"Zh_vs_All" : [(1, 'ZH_vs_All_Zh')],
                 #'BDTsummed_Zh' : [(1, 'BDTsummed_Zh')],
                 #'BDTmultiplied_noZh' : [(1, 'BDTmultiplied_noZh')],
                 #'BDTmultiplied_Zh' : [(1, 'BDTmultiplied_Zh')]
                }
        for discriminant in discriminants.keys() :
            prepareShapes(backgrounds, signals, discriminants[discriminant], discriminant)

def merge_histograms(process, histogram, destination):
    """
    Merge two histograms together. If the destination histogram does not exist, it
    is created by cloning the input histogram

    Parameters:

    process         Name of the current process
    histogram       Pointer to TH1 to merge
    destination     Dict of destination histograms. The key is the current category.

    Return:
    The merged histogram
    """

    if not histogram:
        raise Exception('Missing histogram for %r. This should not happen.' % process)

    #if histogram.GetEntries() == 0:
    #    return

    # Rescale histogram to luminosity, if it's not data
    if process != 'data_obs':
        histogram.Scale(options.luminosity)

    d = destination
    if not d:
        d = histogram.Clone()
        d.SetDirectory(ROOT.nullptr)
    else:
        d = destination
        d.Add(histogram)
    print "Remove negative bins"
    zeroNegativeBins(d)

    return d


def prepareFile(processes_map, categories_map, root_path, discriminant):
    """
    Prepare a ROOT file suitable for Combine Harvester.

    The structure is the following:
      1) Each observable is mapped to a subfolder. The name of the folder is the name of the observable
      2) Inside each folder, there's a bunch of histogram, one per background and signal hypothesis. The name of the histogram is the name of the background.
    """

    import re, hashlib

    print("Preparing ROOT file for combine...")

    output_filename = os.path.join(options.output, 'shapes_MIS_Discri_%s.root' % (discriminant))
    if not os.path.exists(os.path.dirname(output_filename)):
        os.makedirs(os.path.dirname(output_filename))

    hash = hashlib.sha512()
    hash.update(root_path)
    hash.update(output_filename)

    files = [os.path.join(root_path, f) for f in os.listdir(root_path) if f.endswith('_histos.root')]

    # Gather a list of inputs files for each process.
    # The key is the process identifier, the value is a list of files
    # If more than one file exist for a given process, the histograms of each file will
    # be merged together later
    processes_files = {}
    for process, paths in processes_map.items():
        process_files = []
        for path in paths:
            r = re.compile(path, re.IGNORECASE)
            process_files += [f for f in files if r.search(f)]
        processes_files[process] = process_files
        hash.update(process)
        map(hash.update, process_files)

    # Create the list of histograms (nominal + systematics) for each category
    # we are interested in.
    # The key is the category name, and the value is a list of histogram. The list will always
    # contains at least one histogram (the nominal histogram), and possibly more, two per systematic (up & down variation)
    histogram_names = {}
    for category, histogram_name in categories_map.items():
        hash.update(category)
        hash.update(histogram_name)
        r = re.compile(histogram_name, re.IGNORECASE)
        f = ROOT.TFile.Open(processes_files.values()[0][0])
        histogram_names[category] = [n.GetName() for n in f.GetListOfKeys() if r.search(n.GetName())]
        f.Close()

    # Extract list of systematics from the list of histograms derived above
    # This code assumes that *all* categories contains the same systematics (as it should)
    # The systematics list is extracted from the histogram list of the first category
    # The list of expanded histogram name is also extract (ie, regex -> full histogram name)
    systematics = set()
    histograms = {}
    systematics_regex = re.compile('__(.*)(up|down)$', re.IGNORECASE)
    for category, histogram_names in histogram_names.items():
        for histogram_name in histogram_names:
            # FIXME: Untested
            m = systematics_regex.search(histogram_name)
            if m:
                # It's a systematic histogram
                systematics.add(m.group(1))
            else:
                nominal_name = histogram_name
                if category in histograms:
                    # Check that the regex used by the user only math 1 histogram
                    if histograms[category] != nominal_name:
                        raise Exception("The regular expression used for category %r matches more than one histogram: %r and %r" % (category, nominal_name, histograms[category]))
                histograms[category] = nominal_name

    cms_systematics = [CMSNamingConvention(s) for s in systematics]
    for systematic in cms_systematics:
        hash.update(systematic)

    hash = hash.hexdigest()

    #if os.path.exists(output_filename):
    #    # File exists. Check is stored hash is the same as the computed one. If yes, skip the file create
    #    f = ROOT.TFile.Open(output_filename)
    #    stored_hash = f.Get('hash')
    #    if stored_hash and stored_hash.GetTitle() == hash:
    #        print("File %r already exists and contains all the needed shapes. Skipping file generation." % output_filename)
    #        return output_filename, cms_systematics
    #    else:
    #        print("File %r already exists but is no longer up-to-date. It'll be regenerated." % output_filename)

    def dict_get(dict, name):
        if name in dict:
            return dict[name]
        else:
            return None

    # Create final shapes
    shapes = {}
    if options.reweight :
        print "Will reweight some histograms before to feed combine"
    for category, original_histogram_name in histograms.items():
        shapes[category] = {}
        for process, process_files in processes_files.items():
            shapes[category][process] = {}

            for process_file in process_files:
                f = ROOT.TFile.Open(process_file)
                TH1 = f.Get(original_histogram_name)
                if options.reweight :
                    if "DY" in process_file :
                        if not ('DYJetsToLL_M-10to50') in process_file:
                            print "Reweight ", process_file, " by 0.75950" 
                            TH1.Scale(0.75950)
                    if "TTT" in process_file :
                        print "Reweight ", process_file, " by 0.94829" 
                        TH1.Scale(0.94829)
                #shapes[category][process]['nominal'] = merge_histograms(process, f.Get(original_histogram_name), dict_get(shapes[category][process], 'nominal'))
                shapes[category][process]['nominal'] = merge_histograms(process, TH1, dict_get(shapes[category][process], 'nominal'))
                if not process == "data_obs" : 
                    for systematic in systematics:
                        for variation in ['up', 'down']:
                            key = CMSNamingConvention(systematic) + variation.capitalize()
                            shapes[category][process][key] = merge_histograms(process, f.Get(original_histogram_name + '__' + systematic + variation), dict_get(shapes[category][process], key))
                f.Close()

    output_file = ROOT.TFile.Open(output_filename, 'recreate')

    # Store hash
    file_hash = ROOT.TNamed('hash', hash)
    file_hash.Write()

    if options.fake_data:
        for category, processes in shapes.items():
            fake_data = None
            for process, systematics_dict in processes.items():
                #if process.startswith('ggHH'):
                #    continue

                if not fake_data:
                    fake_data = systematics_dict['nominal'].Clone()
                    fake_data.SetDirectory(ROOT.nullptr)
                else:
                    fake_data.Add(systematics_dict['nominal'])

            processes['data_obs'] = {'nominal': fake_data}

    for category, processes in shapes.items():
        output_file.mkdir(category).cd()
        for process, systematics_ in processes.items():
            for systematic, histogram in systematics_.items():
                histogram.SetName(process if systematic == 'nominal' else process + '__' + systematic)
                histogram.Write()
        output_file.cd()

    output_file.Close()
    print("Done. File saved as %r" % output_filename)

    return output_filename, cms_systematics

#def prepareShapesFromBoxes(backgrounds, signals, discriminantName):
#    # Backgrounds is a list of string of the considered backgrounds corresponding to entries in processes_histfactory_mapping
#
#    import CombineHarvester.CombineTools.ch as ch
#
#    if options.file == '' : 
#        print "Must specify input file for box mode"
#        sys.exit()
#    file = options.file 
#
#    categories = [(1, 'yields')]
#
#    cb = ch.CombineHarvester()
#
#    cb.AddObservations(['*'], [''], ['13TeV'], [''], categories)
#
#    cb.AddProcesses(['*'], [''], ['13TeV_2015'], [''], backgrounds, categories, False)
#
#    cb.AddProcesses(['*'], [''], ['13TeV_2015'], [''], signals, categories, True)
#
#    # Systematics
#    #if not options.nosys:
#    #    cb.cp().AddSyst(cb, 'lumi_$ERA', 'lnN', ch.SystMap('era')(['13TeV_2015'], 1.027))
#
#        #cb.cp().AddSyst(cb, '$PROCESS_modeling', 'lnN', ch.SystMap('process')
#        #        (['ttbar'], 1.10)
#        #        (['dy'], 1.30)
#        #        (['SingleTop'], 1.20)
#        #        )
#
#        #cb.cp().AddSyst(cb, '$PROCESS_xsec', 'lnN', ch.SystMap('process')
#        #        (['ttbar'], 1.055842)
#        #        )
#
#    #for systematic in systematics:
#    #    cb.cp().AddSyst(cb, systematic, 'shape', ch.SystMap()(1.00))
#
#    # Import shapes from ROOT file
#    cb.cp().backgrounds().ExtractShapes(file, '$PROCESS_yields', '$PROCESS_yields__$SYSTEMATIC')
#    cb.cp().signals().ExtractShapes(file, '$PROCESS_yields', '$PROCESS_yields__$SYSTEMATIC')
#
#    # Bin by bin uncertainties
#    if options.bbb:
#        print "Making bin by bin uncertainty"
#        bbb = ch.BinByBinFactory()
#        bbb.SetAddThreshold(0.05).SetMergeThreshold(0.5).SetFixNorm(False)
#        bbb.MergeBinErrors(cb.cp().backgrounds())
#        bbb.AddBinByBin(cb.cp().backgrounds(), cb)
#
#    for signal in signals :
#
#        output_prefix = 'MIS_Signal_%s_Discriminant_%s' % (signal, discriminantName)
#
#        output_dir = os.path.join(options.output, '%s' % (signal))
#        if not os.path.exists(output_dir):
#            os.makedirs(output_dir)
#
#        fake_mass = '125'
#
#        # Write card
#        datacard = os.path.join(output_dir, output_prefix + '.dat')
#        cb.cp().mass([fake_mass, "*"]).WriteDatacard(os.path.join(output_dir, output_prefix + '.dat'), os.path.join(output_dir, output_prefix + '_shapes.root'))
#
#        # Write small script to compute the limit
#        workspace_file = os.path.basename(os.path.join(output_dir, output_prefix + '_combine_workspace.root'))
#        script = """#! /bin/bash
#
## If workspace does not exist, create it once
#if [ ! -f {workspace_root} ]; then
#    text2workspace.py {datacard} -m {fake_mass} -o {workspace_root}
#fi
#
## Run limit
#
#combine -H ProfileLikelihood -M Asymptotic -n {name} {workspace_root} -S {systematics}
#""".format(workspace_root=workspace_file, datacard=os.path.basename(datacard), name=output_prefix, fake_mass=fake_mass, systematics=(0 if options.nosys else 1))
#        script_file = os.path.join(output_dir, output_prefix + '_run_limits.sh')
#        with open(script_file, 'w') as f:
#            f.write(script)
#        
#        st = os.stat(script_file)
#        os.chmod(script_file, st.st_mode | stat.S_IEXEC)

def prepareShapes(backgrounds, signals, discriminant, discriminantName):
    # Backgrounds is a list of string of the considered backgrounds corresponding to entries in processes_histfactory_mapping

    import CombineHarvester.CombineTools.ch as ch

    root_path = options.root_path

    categories = discriminant

    def get_hist_regex(r):
        return '^%s(__.*(up|down))?$' % r

    categories_histfactory_mapping = {
            #'BDTsummed_noZh': 'MVA_summed_weight1_withoutZh_kBDT_All_no_cut',
            #'BDTmultiplied_noZh': 'MVA_multiplied_weight1_withoutZh_kBDT_All_no_cut',
            #'TT_weights': 'pp_tt_llbb_tfJetAllEta_minLog_weight_All_hh_llmetjj_HWWleptons_btagM_csv_no_cut',
            #'DY_weights': 'pp_Z_llbb_simple_tfJetAllEta_minLog_weight_All_hh_llmetjj_HWWleptons_btagM_csv_no_cut',
            #'ZH_weights': 'pp_zh_llbb_simple_tfJetAllEta_minLog_weight_All_hh_llmetjj_HWWleptons_btagM_csv_no_cut',
            #'tWp_weights': 'twplus_tfJetAllEta_minLog_weight_All_hh_llmetjj_HWWleptons_btagM_csv_no_cut',
            #'tWm_weights': 'twminus_tfJetAllEta_minLog_weight_All_hh_llmetjj_HWWleptons_btagM_csv_no_cut',
            #'ZZ_weights': 'pp_zz_llbb_simple_tfJetAllEta_minLog_weight_All_hh_llmetjj_HWWleptons_btagM_csv_no_cut',
            #'ATan_DY_TT': 'arcTan_pp_Z_llbb_simple_tfJetAllEta_minus_pp_tt_llbb_tfJetAllEta_All_hh_llmetjj_HWWleptons_btagM_csv_no_cut',
            #'ATan_twp_twm': 'arcTan_twplus_tfJetAllEta_minus_twminus_tfJetAllEta_All_hh_llmetjj_HWWleptons_btagM_csv_no_cut',
            #'csvProd': 'jj_CSVprod_All_hh_llmetjj_HWWleptons_btagM_csv_no_cut',
            #'yields': '^yields(?!(_sf|_df))',
            'yields': get_hist_regex('yields(?!(_sf|_df))'),
            'yields_sf': get_hist_regex('yields_sf'),
            'yields_df': get_hist_regex('yields_df'),
            #'yields': '^yields$',
            #'mll': 'll_M_All_hh_llmetjj_HWWleptons_btagM_csv_no_cut',
            #'mll_ee': 'll_M_ElEl_hh_llmetjj_HWWleptons_btagT_csv_no_cut',
            #'mll_emu': 'll_M_MuEl_hh_llmetjj_HWWleptons_btagT_csv_no_cut',
            #'mll_mumu': 'll_M_MuMu_hh_llmetjj_HWWleptons_btagT_csv_no_cut',
            }
    #if not options.mode == 'box' :
    processes_histfactory_mapping = {
            # Background
            'TT': ['TTTo2L2Nu_13TeV-powheg_Fall15MiniAODv2'],
            #'DY': ['DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_extended_ext0_plus_ext1_plus_ext4', 'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_extended_ext0_plus_ext1_plus_ext3'],
            'DY': ['DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_extended_ext0_plus_ext1_plus_ext4'],
            'DYxx': ['DYxxToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_extended_ext0_plus_ext1_plus_ext4'],
            'DYbx': ['DYbxToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_extended_ext0_plus_ext1_plus_ext4'],
            'DYbb': ['DYbbToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_extended_ext0_plus_ext1_plus_ext4'],
            'DY_10_50': ['DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_extended_ext0_plus_ext1_plus_ext3'],
            #'ttV': ['TT(W|Z)To'],
            'ZH': ['ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8'],
            'ZZ': ['ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_Fall15MiniAODv2'],
            'tWp': ['ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg_Fall15MiniAODv2'],
            'tWm': ['ST_tW_top_5f_inclusiveDecays_13TeV-powheg_Fall15MiniAODv2'],
            'TTV': ['TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Fall15MiniAODv2', 'TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Fall15MiniAODv2', 'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_Fall15MiniAODv2', 'TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_Fall15MiniAODv2'],
            'ST': ['ST_s-channel_4f_leptonDecays_13TeV-amcatnlo_Fall15MiniAODv2', 'ST_t-channel_4f_leptonDecays_13TeV-amcatnlo_Fall15MiniAODv2'],
            'WW' : ['WWTo2L2Nu_13TeV-powheg_Fall15MiniAODv2'],
            'ttH' : ['ttHTobb_M125_13TeV_powheg_pythia8'],
            # Signals
            'TTDM_1_10': ['TTbarDMJets_scalar_Mchi-1_Mphi-10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Fall15MiniAODv2'],
            'TTDM_1_100': ['TTbarDMJets_scalar_Mchi-1_Mphi-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Fall15MiniAODv2'],
            'TTDM_1_500': ['TTbarDMJets_scalar_Mchi-1_Mphi-500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Fall15MiniAODv2'],
            'TTDM_1_10_ps': ['TTbarDMJets_pseudoscalar_Mchi-1_Mphi-10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Fall15MiniAODv2'],
            'TTDM_1_100_ps': ['TTbarDMJets_pseudoscalar_Mchi-1_Mphi-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Fall15MiniAODv2'],
            'TTDM_1_500_ps': ['TTbarDMJets_pseudoscalar_Mchi-1_Mphi-500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Fall15MiniAODv2'],
            
            #'H_200_ZA_50': ['HToZATo2L2B_MH-200_MA-50_13TeV-madgraph_Fall15MiniAODv2'],
            #'H_3000_ZA_2000': ['HToZATo2L2B_MH-3000_MA-2000_13TeV-madgraph_Fall15MiniAODv2'],
            'H_500_ZA_300': ['HToZATo2L2B_MH-500_MA-300_13TeV-madgraph_Fall15MiniAODv2'],
            'H_800_ZA_700': ['HToZATo2L2B_MH-800_MA-700_13TeV-madgraph_Fall15MiniAODv2'],
            
            'Stop_500_325': ['SMS-T2tt_mStop-500_mLSP-325_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Fall15MiniAODv2'],
            'Stop_850_100': ['SMS-T2tt_mStop-850_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Fall15MiniAODv2'],
            'X_400_hh': ['GluGluToRadionToHHTo2B2VTo2L2Nu_M-400_narrow_Fall15MiniAODv2'],
            'X_650_hh': ['GluGluToRadionToHHTo2B2VTo2L2Nu_M-650_narrow_Fall15MiniAODv2'],
            'X_900_hh': ['GluGluToRadionToHHTo2B2VTo2L2Nu_M-900_narrow_Fall15MiniAODv2'],
            'SMhh': ['GluGluToHHTo2B2VTo2L2Nu_node_SM_13TeV-madgraph'],
            #'SMHiggs': ['.*_M125_.*'],
            }
    if not options.fake_data:
        processes_histfactory_mapping['data_obs'] = ['.+Run2015D-16Dec2015.+']
        print "Carefull we use  regex .Run2015D-16Dec2015."
    #else : 
    #    processes_histfactory_mapping = {
    #            # Background
    #            'TT': ['TT_histos.root'],
    #            'DY': ['DY_histos.root'],
    #            #'ttV': ['TT(W|Z)To'],
    #            'ZH': ['ZH_histos.root'],
    #            'VV': ['ZZ_histos.root'],
    #            #'SMHiggs': ['.*_M125_.*'],
    #            }
    #    if not options.fake_data:
    #        processes_histfactory_mapping['data_obs'] = ['data_obs_histos.root']

    file, systematics = prepareFile(processes_histfactory_mapping, categories_histfactory_mapping, root_path, discriminantName)
    
    for signal in signals :
 
        cb = ch.CombineHarvester()

        cb.AddObservations(['*'], [''], ['13TeV'], [''], categories)

        cb.AddProcesses(['*'], [''], ['13TeV_2015'], [''], backgrounds, categories, False)

        cb.AddProcesses(['*'], [''], ['13TeV_2015'], [''], [signal], categories, True)


        # Systematics
        if not options.nosys:
            for systematic in systematics:
                cb.cp().AddSyst(cb, systematic, 'shape', ch.SystMap()(1.00))
            cb.cp().AddSyst(cb, 'lumi_$ERA', 'lnN', ch.SystMap('era')(['13TeV_2015'], 1.027))

            if options.flavSplit :
                print "Bla"
                cb.cp().bin(['yields_df']).AddSyst(cb, '$PROCESS_normalization_df', 'lnN', ch.SystMap('process')
                        (['TT'], 1.001)
                        #(['SingleTop'], 1.)
                        )
                cb.cp().bin(['yields_sf']).AddSyst(cb, '$PROCESS_normalization_sf', 'lnN', ch.SystMap('process')
                        (['TT'], 1.02)
                        (['DY'], 1.04)
                        #(['SingleTop'], 1.)
                        )

            else :
                cb.cp().AddSyst(cb, '$PROCESS_normalization', 'lnN', ch.SystMap('process')
                        (['TT'], 1.02)
                        (['DY'], 1.04)
                        #(['SingleTop'], 1.)
                        )

            #cb.cp().AddSyst(cb, '$PROCESS_xsec', 'lnN', ch.SystMap('process')
            #        (['ttbar'], 1.055842)
            #        )
        if options.nosys and not options.bbb : 
            cb.cp().AddSyst(cb, 'lumi_$ERA', 'lnN', ch.SystMap('era')(['13TeV_2015'], 1.00001)) # Add a negligible systematic (chosen to be lumi) to trick combine

        if options.SF :
            if not options.box : 
                cb.cp().AddSyst(cb, 'SF_$PROCESS', 'rateParam', ch.SystMap('process')
                        (['TT'], 1.)
                        (['DY'], 1.)
                        #(['DYbb'], 1.)
                        #(['DYbx'], 1.)
                        #(['DYxx'], 1.)
                        #(['ZZ'], 1.)
                        #(['ZH'], 1.)
                        #(['tWp'], 1.)
                        #(['tWm'], 1.)
                        )
            else:
                #cb.cp().AddSyst(cb, '$PROCESS_modeling', 'lnN', ch.SystMap('process')
                #        (['TT'], 1.05)
                #        (['DY'], 1.1)
                #        (['ZZ'], 2.0)
                #        (['ZH'], 20.0)
                #        (['tWp'], 2.0)
                #        (['tWm'], 2.0)
                #        )
                cb.cp().AddSyst(cb, 'SF_$PROCESS', 'rateParam', ch.SystMap('process')
                        (['TT'], 1.)
                        (['DY'], 1.)
                        #(['DYbb'], 1.)
                        #(['DYbx'], 1.)
                        #(['DYxx'], 1.)
                        (['ZZ'], 1.)
                        (['ZH'], 1.)
                        (['tWp'], 1.)
                        (['tWm'], 1.)
                        )

        # Import shapes from ROOT file
        cb.cp().backgrounds().ExtractShapes(file, '$BIN/$PROCESS', '$BIN/$PROCESS__$SYSTEMATIC')
        cb.cp().signals().ExtractShapes(file, '$BIN/$PROCESS', '$BIN/$PROCESS__$SYSTEMATIC')

        # Bin by bin uncertainties
        if options.bbb:
            #bbb = ch.BinByBinFactory()
            #bbb.SetAddThreshold(0.05).SetMergeThreshold(0.5).SetFixNorm(False)
            #bbb.MergeBinErrors(cb.cp().backgrounds())
            #bbb.AddBinByBin(cb.cp().backgrounds(), cb)
            bbb = ch.BinByBinFactory()
            if options.SF :
                bbb.SetAddThreshold(0.1).SetMergeThreshold(options.mt).SetFixNorm(False)
                bbb.MergeBinErrors(cb.cp().backgrounds())
            else :
                bbb.SetAddThreshold(0.05).SetMergeThreshold(0.5).SetFixNorm(False)
                #bbb.SetAddThreshold(0.1).SetMergeThreshold(1)
                bbb.MergeBinErrors(cb.cp().backgrounds())
                #bbb_sig = ch.BinByBinFactory()
                #bbb_sig.SetAddThreshold(0.05).SetMergeThreshold(0.5).SetFixNorm(False)
                #bbb_sig.MergeBinErrors(cb.cp().signals())
                #bbb_sig.AddBinByBin(cb.cp().signals(), cb)
                
            bbb.AddBinByBin(cb.cp().backgrounds(), cb)


        output_prefix = 'MIS_Signal_%s_Discriminant_%s' % (signal, discriminantName)

        output_dir = os.path.join(options.output, '%s' % (signal))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        fake_mass = '125'

        # Write card
        datacard = os.path.join(output_dir, output_prefix + '.dat')
        cb.cp().mass([fake_mass, "*"]).WriteDatacard(os.path.join(output_dir, output_prefix + '.dat'), os.path.join(output_dir, output_prefix + '_shapes.root'))

        # Write small script to compute the limit
        workspace_file = os.path.basename(os.path.join(output_dir, output_prefix + '_combine_workspace.root'))
        script = """#! /bin/bash

# If workspace does not exist, create it once
if [ ! -f {workspace_root} ]; then
    text2workspace.py {datacard} -m {fake_mass} -o {workspace_root}
fi

# Run limit

combine -H ProfileLikelihood -M Asymptotic -n {name} {workspace_root} -S {systematics} --rMax 30000
#combine -H ProfileLikelihood -M HybridNew -n {name} {workspace_root} -S {systematics} --testStat LHC --expectedFromGrid 0.5 
#combine -H ProfileLikelihood -M HybridNew -n {name} {workspace_root} -S {systematics} --testStat LHC --expectedFromGrid 0.84 
#combine -H ProfileLikelihood -M HybridNew -n {name} {workspace_root} -S {systematics} --testStat LHC --expectedFromGrid 0.16 
""".format(workspace_root=workspace_file, datacard=os.path.basename(datacard), name=output_prefix, fake_mass=fake_mass, systematics=(0 if options.nosys else 1))
        script_file = os.path.join(output_dir, output_prefix + '_run_limits.sh')
        with open(script_file, 'w') as f:
            f.write(script)
        
        st = os.stat(script_file)
        os.chmod(script_file, st.st_mode | stat.S_IEXEC)

def CMSNamingConvention(syst):
    # Taken from https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWG/HiggsCombinationConventions
    # systlist = ['jec', 'jer', 'elidiso', 'muidiso', 'jjbtag', 'pu', 'trigeff']
    if syst == 'jec':
        return 'CMS_scale_j'
    elif syst == 'jer': 
        return 'CMS_res_j'
    elif syst == 'elidiso': 
        return 'CMS_eff_e'
    elif syst == 'muidiso': 
        return 'CMS_eff_mu'
    elif syst == 'jjbtag': 
        return 'CMS_eff_b'
    elif syst == 'pu': 
        return 'CMS_pu'
    elif syst == 'trigeff': 
        return 'CMS_eff_trigger'
    elif syst == 'pdf':
        return 'pdf'
    elif syst == 'scale':
        return 'QCDscale'
    else:
        return syst
#
# main
#
if __name__ == '__main__':
    main()
