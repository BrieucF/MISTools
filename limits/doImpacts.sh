# /bin/bash
# cd to the directory then 
#call it with doImpacts.sh "combine output file" (do it only with one signal -- mX 400)
text2workspace.py X_400_hh/MIS_Signal_X_400_hh_Discriminant_yields.dat -m 125 -o MIS_Signal_X_400_hh_Discriminant_yields.root
mkdir ImpactsFiles
mkdir Impacts
#pushd ImpactsFiles
combineTool.py -M Impacts -d MIS_Signal_X_400_hh_Discriminant_yields.root -m 125 --robustFit 1 --doInitialFit --rMax 600 -n impacts_combined_asimov_on_data
combineTool.py -M Impacts -d MIS_Signal_X_400_hh_Discriminant_yields.root -m 125 --robustFit 1 --doFits --rMax 600 -n impacts_combined_asimov_on_data --parallel 8
combineTool.py -M Impacts -d MIS_Signal_X_400_hh_Discriminant_yields.root -m 125 -n impacts_combined_asimov_on_data -o Impacts/impacts_combined_asimov_on_data.json --rMax 600 
plotImpacts.py -i Impacts/impacts_combined_asimov_on_data.json -o Impacts/impacts_combined_asimov_on_data
