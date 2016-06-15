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

    if histogram.GetEntries() == 0:
        return

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

    return d

def get_options():
    """
    Parse and return arguments provided by the user
    """
    parser = argparse.ArgumentParser(description='Create shape datacards ready for combine')

    parser.add_argument('-p', '--path', action='store', dest='root_path', type=str, default='/home/fynu/bfrancois/scratch/framework/oct2015/CMSSW_7_4_15/src/cp3_llbb/MISTools/histFactory_mis/mis_8TeV_hist_v3/', help='histFactory output path')

    LUMI = 19700
    parser.add_argument('-l', '--luminosity', action='store', type=float, dest='luminosity', default=LUMI, help='integrated luminosity (default is %f /pb)' % LUMI)

    # hypothesis options
    #parser.add_argument('-c', '--cluster', nargs='+', metavar='CLUSTER', type=str, dest='clusters', default=['SM'], help='clusters hypothesis. Use `all` as an alias for all the clusters')

    # MVA options
    #parser.add_argument('-m', '--trained-cluster', nargs='+', type=str, dest='MVAclusters',
    #    default=['SM'], help='MVA training clusters. Use `all` as an alias for all MVA clusters, or `analysis` to use analysis conditions')

    parser.add_argument('-o', '--output', action='store', dest='output', type=str, default='Cards', help='Output directory')

    parser.add_argument('--fake-data', action='store_true', dest='fake_data',
        help='Use fake data instead of real data')

    parser.add_argument('-s', '--stat', action='store_true', dest='stat_only',
        help='Do not consider systematic uncertainties')

    options = parser.parse_args()
    return options

def main():
    """Main function"""
    global options

    # get the options
    options = get_options()

    backgrounds = ['VV', 'TT', 'DY']
    signals = ['ZH']

    discriminants = { 
            "ThreeBDTTogether_noZh" : [(1, 'TT_vs_All_noZh'), (2, 'DY_vs_All_noZh'), (3, 'ZZ_vs_All_noZh')],
            "ThreeBDTTogether_Zh" : [(1, 'TT_vs_All_Zh'), (2, 'DY_vs_All_Zh'), (3, 'ZZ_vs_All_Zh'), (4, 'ZH_vs_All_Zh')],
            "BDTsummed_noZh" : [(1, 'BDTsummed_noZh')],
            "Zh_vs_All" : [(1, 'ZH_vs_All_Zh')],
            'BDTsummed_Zh' : [(1, 'BDTsummed_Zh')],
            'BDTmultiplied_noZh' : [(1, 'BDTmultiplied_noZh')],
            'BDTmultiplied_Zh' : [(1, 'BDTmultiplied_Zh')]
            }
    for discriminant in discriminants.keys() :
        prepareShapes(backgrounds, signals, discriminants[discriminant], discriminant)

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

    if os.path.exists(output_filename):
        # File exists. Check is stored hash is the same as the computed one. If yes, skip the file create
        f = ROOT.TFile.Open(output_filename)
        stored_hash = f.Get('hash')
        if stored_hash and stored_hash.GetTitle() == hash:
            print("File %r already exists and contains all the needed shapes. Skipping file generation." % output_filename)
            return output_filename, cms_systematics
        else:
            print("File %r already exists but is no longer up-to-date. It'll be regenerated." % output_filename)

    def dict_get(dict, name):
        if name in dict:
            return dict[name]
        else:
            return None

    # Create final shapes
    shapes = {}
    for category, original_histogram_name in histograms.items():
        shapes[category] = {}
        for process, process_files in processes_files.items():
            shapes[category][process] = {}

            for process_file in process_files:
                f = ROOT.TFile.Open(process_file)
                shapes[category][process]['nominal'] = merge_histograms(process, f.Get(original_histogram_name), dict_get(shapes[category][process], 'nominal'))

                for systematic in systematics:
                    for variation in ['up', 'down']:
                        key = CMSNamingConvention(systematic) + variation.capitalize()
                        shapes[category][process][key] = merge_histograms(process, f.Get(original_histogram_name + '__' + systematic + variation), dict_get(shapes[category][process][key]))

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

def prepareShapes(backgrounds, signals, discriminant, discriminantName):
    # Backgrounds is a list of string of the considered backgrounds corresponding to entries in processes_histfactory_mapping

    import CombineHarvester.CombineTools.ch as ch

    root_path = options.root_path

    categories = discriminant

    categories_histfactory_mapping = {
            'BDTsummed_noZh': 'MVA_summed_weight1_withoutZh_kBDT_All_no_cut',
            'BDTmultiplied_noZh': 'MVA_multiplied_weight1_withoutZh_kBDT_All_no_cut',
            'BDTsummed_Zh': 'MVA_summed_weight1_withZh_kBDT_All_no_cut',
            'BDTmultiplied_Zh': 'MVA_multiplied_weight1_withZh_kBDT_All_no_cut',
            'TT_vs_All_noZh': 'MVA_BDT_TT_vs_All_weight1_withoutZh_kBDT_All_no_cut',
            'DY_vs_All_noZh': 'MVA_BDT_DY_vs_All_weight1_withoutZh_kBDT_All_no_cut',
            'ZZ_vs_All_noZh': 'MVA_BDT_ZZ_vs_All_weight1_withoutZh_kBDT_All_no_cut',
            'TT_vs_All_Zh': 'MVA_BDT_TT_vs_All_weight1_withZh_kBDT_All_no_cut',
            'DY_vs_All_Zh': 'MVA_BDT_DY_vs_All_weight1_withZh_kBDT_All_no_cut',
            'ZZ_vs_All_Zh': 'MVA_BDT_ZZ_vs_All_weight1_withZh_kBDT_All_no_cut',
            'ZH_vs_All_Zh': 'MVA_BDT_ZH_vs_All_weight1_withZh_kBDT_All_no_cut',
            }

    processes_histfactory_mapping = {
            # Background
            'TT': ['TT_histos.root'],
            'DY': ['DY_histos.root'],
            #'ttV': ['TT(W|Z)To'],
            'ZH': ['ZH_histos.root'],
            'VV': ['ZZ_histos.root'],
            #'SMHiggs': ['.*_M125_.*'],
            }

    if not options.fake_data:
        processes_histfactory_mapping['data_obs'] = ['data_obs_histos.root']

    file, systematics = prepareFile(processes_histfactory_mapping, categories_histfactory_mapping, root_path, discriminantName)
 
    cb = ch.CombineHarvester()

    cb.AddObservations(['*'], [''], ['8TeV'], [''], categories)

    cb.AddProcesses(['*'], [''], ['8TeV_2015'], [''], backgrounds, categories, False)

    cb.AddProcesses(['*'], [''], ['8TeV_2015'], [''], signals, categories, True)

    # Systematics
    if not options.stat_only:
        cb.cp().AddSyst(cb, 'lumi_$ERA', 'lnN', ch.SystMap('era')(['8TeV_2015'], 1.027))

        #cb.cp().AddSyst(cb, '$PROCESS_modeling', 'lnN', ch.SystMap('process')
        #        (['ttbar'], 1.10)
        #        (['dy'], 1.30)
        #        (['SingleTop'], 1.20)
        #        )

        #cb.cp().AddSyst(cb, '$PROCESS_xsec', 'lnN', ch.SystMap('process')
        #        (['ttbar'], 1.055842)
        #        )

    for systematic in systematics:
        cb.cp().AddSyst(cb, systematic, 'shape', ch.SystMap()(1.00))

    # Import shapes from ROOT file
    cb.cp().backgrounds().ExtractShapes(file, '$BIN/$PROCESS', '$BIN/$PROCESS__$SYSTEMATIC')
    cb.cp().signals().ExtractShapes(file, '$BIN/$PROCESS', '$BIN/$PROCESS__$SYSTEMATIC')

    # Bin by bin uncertainties
    if not options.stat_only:
        bbb = ch.BinByBinFactory()
        bbb.SetAddThreshold(0.05).SetMergeThreshold(0.5).SetFixNorm(False)
        bbb.MergeBinErrors(cb.cp().backgrounds())
        bbb.AddBinByBin(cb.cp().backgrounds(), cb)

    for signal in signals :

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

combine -H ProfileLikelihood -M Asymptotic -n {name} {workspace_root} -S {systematics}
""".format(workspace_root=workspace_file, datacard=os.path.basename(datacard), name=output_prefix, fake_mass=fake_mass, systematics=(0 if options.stat_only else 1))
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
