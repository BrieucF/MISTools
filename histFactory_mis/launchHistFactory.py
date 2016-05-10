#!/nfs/soft/python/python-2.7.5-sl6_amd64_gcc44/bin/python

# usage dans histfactory_hh : python launchHistFactory.py outputName [submit]

import sys, os, json
sys.path.append("../../CommonTools/histFactory/")
import copy
import datetime

import argparse

parser = argparse.ArgumentParser(description='Facility to submit histFactory jobs on condor.')
parser.add_argument('-o', '--output', dest='output', default=str(datetime.date.today()), help='Name of the output directory.')
parser.add_argument('-p', '--plotter', dest='plotter', default="generatePlots.py", help='Code generating the plots.')
parser.add_argument('-r', '--remove', help='Overwrite output directory if it already exists.', action="store_true")
parser.add_argument('--skip', help='Skip the building part.', action="store_true")
parser.add_argument('--noExec', help='Skip the execution part.', action="store_true")
parser.add_argument('-j', '--json', dest='json', default='/home/fynu/bfrancois/scratch/llbb_release/modifieForJetEtaTooHigh/CMSSW_5_3_20/src/UserCode/ttbar_effth_delphes/analyzer/jsons/mis_8TeV.json', help='Directory of the json with input files.')

args = parser.parse_args()

jsonName = args.json 
jsonFile = open(jsonName)
datasetDict = json.load(jsonFile)
for datasetName in datasetDict.keys() :
    # Carefull, must not be data!
    rootFileName = datasetDict[datasetName]["files"][0]
    break

if os.path.isdir(args.output) :
    if args.remove :
        print "Are you sure you want to execute the following command ?"
        print "rm -r " + args.output
        print "Type enter if yes, ctrl-c if not."
        raw_input()
        os.system("rm -r " + args.output)
        print "Deleted ", args.output, " folder."
    elif not args.skip :
        print "Output directory : %s exists, exiting..."%args.output
        sys.exit()

print rootFileName
if not args.skip :
    os.system("../../CommonTools/histFactory/build/createPlotter.sh %s %s %s"%(rootFileName, args.plotter, args.output))
if not args.noExec :
    currentPath = os.getcwd()
    os.chdir(args.output)
    os.system("./build/parallelizedPlotter.py -j 4 %s"%(args.json))
    os.chdir(currentPath)
