== general ===
The photonproduction of BH crossection is based on formula Eur.Phys.J.C23:675-689,2002
Then virtual photon flux is calculated by equivalent photon approximation (EPA) by using formula 3 of
http://lss.fnal.gov/archive/other/lpc-94-35.pdf 

anaTCS.C can be used to analyze the output tree
It's too slow to run as script, better run as binary
root 'anaTCS.C+("../tcs_gen.root","SoLID",false)'
choose output root file
choose SoLID or CLAS12 detector
use "true" or "false" to turn on or off momentum smearing 
The total event counts is shown on screen

TCDGen only simulate photon 4-momentum going along z direction and its Q2 cut only affect photon counts
momentum smearing are just based some detector parameters and we can do further Q2 cut

root file with acceptance histogram is needed for anaTCS.C to work
for SoLID JPsi setup, "acceptance_solid_JPsi_electron_target315_output.root"
for CLAS12 setup, "clasev_acceptance.root"

== case study =====

Comparing to ~3000 event counts from CLAS12 RGA Fall 2008 data used in Pierre's thesis 

Use input file "GenOptions_clas12.dat" and cut in "anaTCS.C" to ensure the following condition
"Eb=10.6,Q2<0.15,1.5<Q'<3,4<Eg<10.6,0.15<-t<0.8,integrated lumi 200/fb,geometry acceptance by clas12 fast simulation,
overall detection eff 50%"
running "root 'anaTCS.C+("../tcs_gen.root","CLAS12",false)'" will get total event counts 19034 which is 6 times of 3000
It could be from inaccurate geometry acceptance and or overall detection eff?

