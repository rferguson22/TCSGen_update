The photonproduction of BH crossection is based on formula Eur.Phys.J.C23:675-689,2002
Then virtual photon flux is calculated by equivalent photon approximation (EPA) by using formula 3 of
http://lss.fnal.gov/archive/other/lpc-94-35.pdf 

anaTCS.C can be used to analyze the output tree
It's too slow to run as script, better run as binary
root 'anaTCS.C+("../tcs_gen.root","SoLID",false)'
choose output root file
choose SoLID or CLAS12 detector
use "true" or "false" to turn on or off momentum smearing

root file with acceptance histogram is needed for anaTCS.C to work
for SoLID JPsi, "wget https://github.com/JeffersonLab/solid_gemc/raw/master/analysis/acceptance/result_JPsi/201501/acceptance_solid_JPsi_electron_target315_output.root"

momentum smearing are just based some parameters in code