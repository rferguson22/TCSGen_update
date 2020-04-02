// #include <sys/time.h>
// #include <sys/resource.h>
// #include <vector>
#include <iostream>
#include <TMath.h> 
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TRandom.h>
#include <TApplication.h>

using namespace std;

int anaTCS(string input_filename,string detector,bool Is_smear){
// int main(int argc,char *argv[]){
//  string  input_filename=argv[1];
//   string  detector=argv[2];

//     const rlim_t kStackSize = 500L * 1024L * 1024L;   // min stack size = 64 Mb
//     struct rlimit rl;
//     int result;
// 
//     result = getrlimit(RLIMIT_STACK, &rl);
//     if (result == 0)
//     {
// 	cout << "rlim_cur " <<  rl.rlim_cur << endl;
//         if (rl.rlim_cur < kStackSize)
//         {
//             rl.rlim_cur = kStackSize;
//             result = setrlimit(RLIMIT_STACK, &rl);
//             if (result != 0)
//             {
//                 fprintf(stderr, "setrlimit returned result = %d\n", result);
//             }
//         }
// 	cout << "rlim_cur " <<  rl.rlim_cur << endl;        
//     }
  
      if(Is_smear) cout << "do smearing according to which detector" << endl;
      else cout << "do NO smearing" << endl;

// gROOT->Reset();
gStyle->SetPalette(1);
// gStyle->SetOptStat(1111111);
gStyle->SetOptStat(0);


const double DEG=180/3.1416;

//   Double_t cov= 1e-12 * 1e-24; //pb to cm2 coversion
//   Double_t lumi = 1e37;  // 1e37/cm2/s is from 3nA on 15cm long LH2 target
//   Double_t br = 1.;
//   Double_t range_angle = 4*3.1415926;
//   Double_t range_P= (3.62-0)*(1.6-0);  // Q' and t range
//   Double_t time = 50*3600*24;  //50 days in seconds
//   Double_t eff = 0.85;
//   Double_t overall_NOneve = cov *lumi * br * range_angle * range_P * eff * time;  

TLorentzVector target;
double Pmax,smin,smax,etamin,etamax,tmin,tmax,thetamin,thetamax;
double resolution_p[2],resolution_theta[2],resolution_phi[2];  // 0 is forward, 1 is large angle

int index;

TH2F *hacceptance_negative_largeangle,*hacceptance_negative_forwardangle,*hacceptance_positive_largeangle,*hacceptance_positive_forwardangle,*hacceptance_negative,*hacceptance_positive;
TH3F *hacceptance_PThetaPhi_negative,*hacceptance_PThetaPhi_positive;
if (detector=="SoLID"){  
// TFile *acceptancefile_negative=new TFile("acceptance_solid_CLEO_JPsi_negative_output.root");
// TH2F *hacceptance_6_negative=(TH2F*) acceptancefile_negative->Get("acceptance_6");
// TH2F *hacceptance_7_negative=(TH2F*) acceptancefile_negative->Get("acceptance_7");
// TFile *acceptancefile_positive=new TFile("acceptance_solid_CLEO_JPsi_positive_output.root");
// TH2F *hacceptance_6_positive=(TH2F*) acceptancefile_positive->Get("acceptance_6");
// TH2F *hacceptance_7_positive=(TH2F*) acceptancefile_positive->Get("acceptance_7");

// TFile *acceptancefile_negative=new TFile("acceptance_solid_CLEO_JPsi_negative_target300_output.root");
// hacceptance_negative_largeangle=(TH2F*) acceptancefile_negative->Get("acceptance_largeangle");
// hacceptance_negative_forwardangle=(TH2F*) acceptancefile_negative->Get("acceptance_forwardangle");
// hacceptance_positive_largeangle=(TH2F*) acceptancefile_negative->Get("acceptance_largeangle");
// hacceptance_positive_forwardangle=(TH2F*) acceptancefile_negative->Get("acceptance_forwardangle");
  
TFile *acceptancefile_negative=new TFile("acceptance_solid_JPsi_electron_target315_output.root");  
hacceptance_negative_largeangle=(TH2F*) acceptancefile_negative->Get("acceptance_ThetaP_largeangle");
hacceptance_negative_forwardangle=(TH2F*) acceptancefile_negative->Get("acceptance_ThetaP_forwardangle");
hacceptance_positive_largeangle=(TH2F*) acceptancefile_negative->Get("acceptance_ThetaP_largeangle");
hacceptance_positive_forwardangle=(TH2F*) acceptancefile_negative->Get("acceptance_ThetaP_forwardangle");
 
// TH2F *hacceptance_negative;
// hacceptance_negative->Add(hacceptance_negative_largeangle,hacceptance_negative_forwardangle);
// TH2F *hacceptance_positive;
// hacceptance_positive->Add(hacceptance_positive_largeangle,hacceptance_positive_forwardangle);

hacceptance_negative=(TH2F*) hacceptance_negative_largeangle->Clone();
hacceptance_negative->Add(hacceptance_negative_forwardangle);
hacceptance_positive=(TH2F*) hacceptance_positive_largeangle->Clone();
hacceptance_positive->Add(hacceptance_positive_forwardangle);

//   TCanvas *c_acceptance = new TCanvas("acceptance","acceptance",1200,900);
//   c_acceptance->Divide(2,3);
//   c_acceptance->cd(1);
//   hacceptance_negative_largeangle->Draw("colz");
//   c_acceptance->cd(3);
//   hacceptance_negative_forwardangle->Draw("colz");  
//   c_acceptance->cd(5);
//   hacceptance_negative->Draw("colz");   
//   c_acceptance->cd(2);
//   hacceptance_positive_largeangle->Draw("colz");
//   c_acceptance->cd(4);
//   hacceptance_positive_forwardangle->Draw("colz");
//   c_acceptance->cd(6);
//   hacceptance_positive->Draw("colz");

//   TCanvas *c_acceptance_ratio = new TCanvas("acceptance_ratio","acceptance_ratio",1000,900);  
//   TH2F *hacceptance_ratio=(TH2F*) hacceptance_negative->Clone();
//   hacceptance_ratio->Divide(hacceptance_negative,hacceptance_positive);
//   hacceptance_ratio->Draw("colz"); 
  
  target.SetPxPyPzE(0.,0.,0.,0.938);
  Pmax=3;
  smin=10;
  smax=25;
  etamin=0.;
  etamax=0.5;
  tmin=0;  
  tmax=4;
  index=4;
  thetamin=0;  
  thetamax=40;
//resolution mom 2%, theta 0.6mr, phi 5mr   
  resolution_p[0]=0.03;resolution_theta[0]=1e-3;resolution_phi[0]=6e-3;
  resolution_p[1]=0.02;resolution_theta[1]=0.6e-3;resolution_phi[1]=5e-3;   
}
else if (detector=="CLAS12"){
//   TFile *acceptancefile=new TFile("clas12_acceptance_pipele.root");
//   TFile *acceptancefile=new TFile("clas12_acceptance_pipele_largebin.root");
  TFile *acceptancefile=new TFile("clasev_acceptance.root");
  
  hacceptance_PThetaPhi_positive=(TH3F*) acceptancefile->Get("acceptance_PThetaPhi_pip");  
  hacceptance_PThetaPhi_negative=(TH3F*) acceptancefile->Get("acceptance_PThetaPhi_ele");
  
//   TCanvas *c_acceptance = new TCanvas("acceptance","acceptance",1200,900);
//   c_acceptance->Divide(2,1);
//   c_acceptance->cd(1);
//   hacceptance_PThetaPhi_positive->Draw();
//   c_acceptance->cd(2);
//   hacceptance_PThetaPhi_negative->Draw();  
  
  target.SetPxPyPzE(0.,0.,0.,0.938);  
  Pmax=3; 
  smin=10;
  smax=25;
  etamin=0.;  
  etamax=0.5;
  tmin=0;    
  tmax=4;
  index=4; 
  thetamin=0;    
  thetamax=60;  
//resolution mom 1%, theta 1mr, phi 4mr   
  resolution_p[0]=0.01;resolution_theta[0]=1e-3;resolution_phi[0]=4e-3;
  resolution_p[1]=0.05;resolution_theta[1]=10e-3;resolution_phi[1]=5e-3;   
}
else if (detector=="EIC"){
  target.SetPxPyPzE(0., -60*sin(6./DEG), -60*cos(6./DEG), sqrt(60*60+0.938*0.938));
  Pmax=70;
  smin=0;
  smax=4000;
  etamin=0.;
  etamax=1;
  tmin=0;  
  tmax=300;
  index=1;
  thetamin=0;    
  thetamax=180;  
//resolution mom 4%, theta 1mr, phi 6mr   
  resolution_p[0]=0.01;resolution_theta[0]=1e-3;resolution_phi[0]=4e-3;
  resolution_p[1]=0.05;resolution_theta[1]=10e-3;resolution_phi[1]=5e-3;   
}
else {
  cout << " I don't know your detector " << endl;
  return 0;
}

char output_filename[200];
// sprintf(output_filename, "%s_%s_output.root",input_filename.substr(input_filename.find("/"),input_filename.rfind(".")).c_str(),detector.c_str());
sprintf(output_filename, "output.root");
TFile *outputfile=new TFile(output_filename, "recreate");

TH1F *hcrs_BH_log=new TH1F("crs_BH_log","crs_BH_log",24,-12,12);

TH1F *hcrs_t_stat=new TH1F("crs_t_stat","crs_t_stat",5,0.1,1.1);
TH1F *hcrs_t=new TH1F("crs_t","crs_t",5,0.1,1.1);
// TH1F *hcrs_t_paper=new TH1F("crs_t_paper","crs_t_paper",40,0.,0.4);
// TH1F *hcrs_t_paper_nopsf=new TH1F("crs_t_paper_nopsf","crs_t_paper_nopsf",40,0.,0.4);
TH1F *hcrs_t_paper=new TH1F("crs_t_paper","crs_t_paper",100,0.,1);
TH1F *hcrs_t_paper_nopsf=new TH1F("crs_t_paper_nopsf","crs_t_paper_nopsf",100,0.,1);



TH1F *hEg=new TH1F("Eg","Eg",120,0,12);
TH1F *hflux_factor=new TH1F("flux_factor","flux_factor",200,-0.001,0.003);
TH1F *hEgflux_factor=new TH1F("Egflux_factor","Egflux_factor",120,0,12);
TH2F *hEg_flux_factor=new TH2F("Eg_flux_factor","Eg_flux_factor",120,0,12,200,-0.001,0.003);

const int n=5;

TH2F *heta_xi[n];
TH3F *ht_Q2_xi[n];
TH2F *hQ2_xi[n];
TH1F *hphoton_m2[n];
TH2F *hphoton_theta_mom[n],*hproton_theta_mom[n],*helectron_theta_mom[n],*hpositron_theta_mom[n];
TH2F *helectron_positron_theta_mom_ratio[n];
TH2F *hep_mom[n],*hep_theta[n];
TH1F *ht[n],*hs[n],*hep_InvM[n];
TH2F *hQ2_t[n];
TH1F *htau[n];
TH1F *heta[n];
TH2F *htheta_phi_CM[n];
TH1F *hphi_CM[n];
TH2F *htheta_phi_CM_bin[4][n];
TH1F *hphi_CM_bin[4][n];
TH2F *hMissP_MM2[n],*hQ2_MM2[n],*hMissPxPy[n];
TH1F *hMissMM2[n],*hMissPt[n];
for(int k=0;k<n;k++){  
  char hstname[100];

  sprintf(hstname,"eta_xi_%i",k);
  heta_xi[k]=new TH2F(hstname,hstname,100,0,1,100,0,1);
  sprintf(hstname,"t_Q2_xi_%i",k);
  ht_Q2_xi[k]=new TH3F(hstname,hstname,100,0,0.5,70,3,10,100,tmin,tmax);
  ht_Q2_xi[k]->SetTitle(";#xi;Q'2(GeV^{2});t(GeV^{2})");      
  sprintf(hstname,"Q2_xi_%i",k);
  hQ2_xi[k]=new TH2F(hstname,hstname,100,0,0.5,70,3,10);
  hQ2_xi[k]->SetTitle(";#xi;Q'2(GeV^{2})");
  sprintf(hstname,"Q2_t_%i",k);  
  hQ2_t[k]=new TH2F(hstname,hstname,100,tmin,tmax,70,3,10);
  hQ2_t[k]->SetTitle(";t(GeV^{2});Q'2(GeV^{2})");
    
  sprintf(hstname,"photon_m2_%i",k);
  hphoton_m2[k]=new TH1F(hstname,hstname,100,-1,1);

  sprintf(hstname,"photon_theta_mom_%i",k);
  hphoton_theta_mom[k]=new TH2F(hstname,hstname,30,0,3,110,0,11);
  hphoton_theta_mom[k]->SetTitle("Photon;#theta (degree);Mom (GeV)");    
  sprintf(hstname,"proton_theta_mom_%i",k);
  hproton_theta_mom[k]=new TH2F(hstname,hstname,180,0,180,110,0,Pmax);
  hproton_theta_mom[k]->SetTitle("Proton;#theta (degree);Mom (GeV)");  
  sprintf(hstname,"electron_theta_mom_%i",k);
  helectron_theta_mom[k]=new TH2F(hstname,hstname,180,0,180,110,0,11);  
  helectron_theta_mom[k]->SetTitle("electron;#theta (degree);Mom (GeV)");  
  sprintf(hstname,"positron_theta_mom_%i",k);
  hpositron_theta_mom[k]=new TH2F(hstname,hstname,180,0,180,110,0,11);
  hpositron_theta_mom[k]->SetTitle("positron;#theta (degree);Mom (GeV)");  
  
  sprintf(hstname,"electron_positron_theta_mom_ratio_%i",k);  
  helectron_positron_theta_mom_ratio[k]=new TH2F(hstname,hstname,180,0,180,110,0,11);
    
  sprintf(hstname,"ep_mom_%i",k);
  hep_mom[k]=new TH2F(hstname,hstname,55,0,11,55,0,11);
  hep_mom[k]->SetTitle("electron Mom vs positron Mom;Mom (GeV);Mom (GeV)");    
  sprintf(hstname,"ep_theta_%i",k);
  hep_theta[k]=new TH2F(hstname,hstname,180,0,180,180,0,180);
  hep_theta[k]->SetTitle("electron #theta vs positron #theta;#theta (degree);#theta (degree)");        
  sprintf(hstname,"t_%i",k);
  ht[k]=new TH1F(hstname,hstname,50,tmin,tmax);
  sprintf(hstname,"s_%i",k);
  hs[k]=new TH1F(hstname,hstname,50,smin,smax);
  sprintf(hstname,"ep_InvM_%i",k);
  hep_InvM[k]=new TH1F(hstname,hstname,100,0,4);
  
  sprintf(hstname,"tau_%i",k);
  htau[k]=new TH1F(hstname,hstname,100,0,1);
  sprintf(hstname,"eta_%i",k);
  heta[k]=new TH1F(hstname,hstname,100,0,1);
  
  sprintf(hstname,"theta_phi_CM_%i",k);
  htheta_phi_CM[k]=new TH2F(hstname,hstname,180,-3.15,3.15,180,0,3.15);
  sprintf(hstname,"phi_CM_%i",k);
  hphi_CM[k]=new TH1F(hstname,hstname,180,-3.15,3.15);
  hphi_CM[k]->SetTitle(";#phi_CM;d#sigma/(dQ dt d#theta d#phi) (pb/GeV^3)");  

  sprintf(hstname,"MissP_MM2_%i",k);
  hMissP_MM2[k]=new TH2F(hstname,hstname,100,-0.1,0.1,500,0,5);
  hMissP_MM2[k]->SetTitle(";MM^{2} (GeV^{2});Missing P(GeV)");    
  sprintf(hstname,"Q2_MM2_%i",k);
  hQ2_MM2[k]=new TH2F(hstname,hstname,100,-0.1,0.1,100,-1,1);
  hQ2_MM2[k]->SetTitle(";MM^{2} (GeV^{2});Q^{2} (GeV^{2})");      
  sprintf(hstname,"MissPxPy_%i",k);
  hMissPxPy[k]=new TH2F(hstname,hstname,100,-1,1,100,-1,1);
  hMissPxPy[k]->SetTitle(";Missing P_{x} (GeV);Missing P_{y} (GeV)");      
  
  sprintf(hstname,"MissMM2_%i",k);
  hMissMM2[k]=new TH1F(hstname,hstname,100,-0.1,0.1);
  sprintf(hstname,"hMissPt_%i",k);
  hMissPt[k]=new TH1F(hstname,hstname,100,1,1);
  
  for(int l=0;l<4;l++){  
    sprintf(hstname,"theta_phi_CM_bin_%i_%i",l,k);
    htheta_phi_CM_bin[l][k]=new TH2F(hstname,hstname,180,-3.15,3.15,180,0,3.15);
    sprintf(hstname,"phi_CM_bin_%i_%i",l,k);    
    hphi_CM_bin[l][k]=new TH1F(hstname,hstname,180,-3.15,3.15);    
    hphi_CM_bin[l][k]->SetTitle(";#phi_CM;d#sigma/(dQ dt d#theta d#phi) (pb/GeV^3)");
  }
  
}

const int Nsbin=6;
const int Netabin=6;
const int NQ2bin=3;
const int Ntbin=4;
double etabin_edge[Netabin+1]={0.100,0.140,0.175,0.195,0.210,0.230,0.400};
double sbin_edge[Nsbin+1]={11,14.5,16,17,18,19.5,23};
double Q2bin_edge[NQ2bin+1]={4.0,4.4,4.9,9.0};
// double Q2bin_edge[NQ2bin+1]={4,4.25,4.75,5.25,6.0,9};
// double tbin_edge[Ntbin+1]={0.0,0.2,0.4,0.8,1.6,3.2};
double tbin_edge[Ntbin+1]={0.0,0.1,0.3,0.7,1.5};
// double tbin_edge[Ntbin+1]={0.0,0.3,0.34,0.38,0.42,0.8,1.6,3.2};

const int NQ2bin_etabin=5;
// double Q2bin_etabin_edge[Netabin][NQ2bin_etabin+1]={
double Q2bin_etabin_edge[36]={
4.00,4.15,4.25,4.40,4.60,5.00,
4.00,4.25,4.45,4.70,5.05,6.00,
4.00,4.10,4.35,4.70,5.40,6.50,
4.00,4.20,4.40,4.65,5.25,7.00,
4.00,4.35,4.65,4.85,5.25,7.50,
4.00,4.65,5.20,5.60,6.00,9.00
};

cout << "Q2bin_etabin_center" <<  endl;
for(int etabin=0;etabin<Netabin;etabin++){
for(int Q2bin=0;Q2bin<NQ2bin_etabin;Q2bin++){
  cout << (Q2bin_etabin_edge[etabin*(NQ2bin_etabin+1)+Q2bin+1]+Q2bin_etabin_edge[etabin*(NQ2bin_etabin+1)+Q2bin])/2 << "\t";
}
  cout << endl;
}

double countbin_etabin[Netabin][NQ2bin_etabin][1];
for(int etabin=0;etabin<Netabin;etabin++){
for(int Q2bin=0;Q2bin<NQ2bin_etabin;Q2bin++){
for(int tbin=0;tbin<1;tbin++){ 
  countbin_etabin[etabin][Q2bin][tbin]=0;
}}}

TH2F *htheta_phi_CM_etabin[Netabin];
TH2F *ht_Q2_etabin[Netabin],*ht_Q2_etabin_even[Netabin];
for(int etabin=0;etabin<Netabin;etabin++){
    char hstname[100];
  sprintf(hstname,"t_Q2_etabin_%i",etabin);  
  ht_Q2_etabin[etabin]=new TH2F(hstname,hstname,120,0,1,140,3,10);     
  sprintf(hstname,"t_Q2_etabin_even_%i",etabin);  
  ht_Q2_etabin_even[etabin]=new TH2F(hstname,hstname,120,0,1,140,3,10);     
  sprintf(hstname,"theta_phi_CM_etabin_%i",etabin);  
  htheta_phi_CM_etabin[etabin]=new TH2F(hstname,hstname,180,-3.15,3.15,180,0,3.15);
}


const int Netabin_Q2bin=12;
double etabin_Q2bin_edge[39]={
0.100,0.125,0.135,0.145,0.155,0.175,0.180,0.190,0.195,0.202,0.207,0.220,0.280,
0.115,0.135,0.140,0.150,0.158,0.170,0.180,0.195,0.202,0.210,0.220,0.230,0.280,
0.125,0.155,0.170,0.180,0.190,0.200,0.215,0.225,0.232,0.240,0.250,0.265,0.380
};

cout << "etabin_Q2bin_center" <<  endl;
for(int Q2bin=0;Q2bin<NQ2bin;Q2bin++){
for(int etabin=0;etabin<Netabin_Q2bin;etabin++){
  cout << (etabin_Q2bin_edge[Q2bin*(Netabin_Q2bin+1)+etabin+1]+etabin_Q2bin_edge[Q2bin*(Netabin_Q2bin+1)+etabin])/2 << "\t";
}
  cout << endl;
}


double countbin_Q2bin[NQ2bin][Netabin_Q2bin][1];
for(int Q2bin=0;Q2bin<NQ2bin;Q2bin++){
for(int etabin=0;etabin<Netabin_Q2bin;etabin++){
for(int tbin=0;tbin<1;tbin++){ 
  countbin_Q2bin[Q2bin][etabin][tbin]=0;
}}}

TH2F *htheta_phi_CM_Q2bin[NQ2bin];
TH2F *ht_eta_Q2bin[NQ2bin],*ht_eta_Q2bin_even[NQ2bin];
for(int Q2bin=0;Q2bin<NQ2bin;Q2bin++){
    char hstname[100];
  sprintf(hstname,"t_eta_Q2bin_%i",Q2bin);  
  ht_eta_Q2bin[Q2bin]=new TH2F(hstname,hstname,120,0,1,100,0,0.5);     
  sprintf(hstname,"t_eta_Q2bin_even_%i",Q2bin);  
  ht_eta_Q2bin_even[Q2bin]=new TH2F(hstname,hstname,120,0,1,100,0,0.5);     
  sprintf(hstname,"theta_phi_CM_Q2bin_%i",Q2bin);  
  htheta_phi_CM_Q2bin[Q2bin]=new TH2F(hstname,hstname,180,-3.15,3.15,180,0,3.15);
}

TH2F *ht_Q2_sbin[Nsbin],*ht_Q2_sbin_even[Nsbin];
for(int sbin=0;sbin<Nsbin;sbin++){
  char hstname[100];
  sprintf(hstname,"t_Q2_sbin_%i",sbin);  
  ht_Q2_sbin[sbin]=new TH2F(hstname,hstname,120,0,1,140,3,10);     
  sprintf(hstname,"t_Q2_sbin_even_%i",sbin);  
  ht_Q2_sbin_even[sbin]=new TH2F(hstname,hstname,120,0,1,140,3,10);    
}


// double countbin_etabin[Netabin][Ntbin][NQ2bin];
// // TH2F *hThetaPhiCM_etabin_Q2bin_tbin[Netabin][NQ2bin][Ntbin]; 
// for(int etabin=0;etabin<Netabin;etabin++){
// for(int Q2bin=0;Q2bin<NQ2bin;Q2bin++){
// for(int tbin=0;tbin<Ntbin;tbin++){  
// //   char hstname[100];
// //   sprintf(hstname,"ThetaPhiCM_etabin_Q2bin_tbin_%i_%i_%i",etabin,Q2bin,tbin);  
// //   hThetaPhiCM_etabin_Q2bin_tbin[etabin][Q2bin][tbin]=new TH2F(hstname,hstname,45,-3.15,3.15,45,0,3.15);  
//   countbin_etabin[etabin][Q2bin][tbin]=0;
// }}}

// vector<TH2F> *hThetaPhiCM_etabin_Q2bin_tbin;
// for(int etabin=0;etabin<Netabin;etabin++){
// for(int tbin=0;tbin<Ntbin;tbin++){
// for(int Q2bin=0;Q2bin<NQ2bin;Q2bin++){
//   char hstname[100];
//   sprintf(hstname,"ThetaPhiCM_etabin_Q2bin_tbin_%i_%i_%i",etabin,tbin,Q2bin);  
//   TH2F *hThetaPhiCM_etabin_Q2bin_tbin_temp=new TH2F(hstname,hstname,90,-3.15,3.15,45,0,3.15);  
//   hThetaPhiCM_etabin_Q2bin_tbin->push_back(*hThetaPhiCM_etabin_Q2bin_tbin_temp);
// }}}

TH1F *hs_final=new TH1F("s_final","s_final",50,smin,smax);
TH1F *ht_final=new TH1F("t_final","t_final",50,tmin,tmax);
TH1F *heta_final=new TH1F("eta_final","eta_final",50,etamin,etamax);
TH1F *hQ2_final=new TH1F("Q2_final","Q2_final",70,3,10);

TH2F *hs_t_final=new TH2F("s_t_final","s_t_final",50,smin,smax,50,tmin,tmax);
TH2F *hs_Q2_final=new TH2F("s_Q2_final","s_Q2_final",50,smin,smax,70,3,10);
TH2F *heta_t_final=new TH2F("eta_t_final","eta_t_final",50,etamin,etamax,50,tmin,tmax);
TH2F *heta_Q2_final=new TH2F("eta_Q2_final","eta_Q2_final",50,etamin,etamax,70,3,10);
TH2F *ht_Q2_final=new TH2F("t_Q2_final","t_Q2_final",50,tmin,tmax,70,3,10);

TH3F *heta_t_Q2_final=new TH3F("eta_t_Q2_final","eta_t_Q2_final",50,etamin,etamax,50,tmin,tmax,70,3,10);
TH3F *hs_t_Q2_final=new TH3F("s_t_Q2_final","s_t_Q2_final",50,smin,smax,50,tmin,tmax,70,3,10);

TH1F *hs_final_bin=new TH1F("s_final_bin","s_final_bin",Nsbin,sbin_edge);
TH1F *ht_final_bin=new TH1F("t_final_bin","t_final_bin",Ntbin,tbin_edge);
TH1F *heta_final_bin=new TH1F("eta_final_bin","eta_final_bin",Netabin,etabin_edge);
TH1F *hQ2_final_bin=new TH1F("Q2_final_bin","Q2_final_bin",NQ2bin,Q2bin_edge);


TH2F *hs_t_final_bin=new TH2F("s_t_final_bin","s_t_final_bin",Nsbin,sbin_edge,Ntbin,tbin_edge);
TH2F *hs_Q2_final_bin=new TH2F("s_Q2_final_bin","s_Q2_final_bin",Nsbin,sbin_edge,NQ2bin,Q2bin_edge);
TH2F *heta_t_final_bin=new TH2F("eta_t_final_bin","eta_t_final_bin",Netabin,etabin_edge,Ntbin,tbin_edge);
TH2F *heta_Q2_final_bin=new TH2F("eta_Q2_final_bin","eta_Q2_final_bin",Netabin,etabin_edge,NQ2bin,Q2bin_edge);
TH2F *ht_Q2_final_bin=new TH2F("t_Q2_final_bin","t_Q2_final_bin",Ntbin,tbin_edge,NQ2bin,Q2bin_edge);

TH3F *heta_t_Q2_final_bin=new TH3F("eta_t_Q2_final_bin","eta_t_Q2_final_bin",Netabin,etabin_edge,Ntbin,tbin_edge,Ntbin,tbin_edge);
TH3F *hs_t_Q2_final_bin=new TH3F("s_t_Q2_final_bin","s_t_Q2_final_bin",Nsbin,sbin_edge,Ntbin,tbin_edge,Ntbin,tbin_edge);

TH2F *hgen_t_Q2[20];
TH2F *hgen_theta_phi_CM[20];
for(int k=0;k<20;k++){
   char hstname[100];
//    sprintf(hstname,"gen_t_InvM2_ep_%i",k);
//    hgen_t_Q2[k]=new TH2F(hstname,hstname,200,0,20,900,3,12);
      sprintf(hstname,"gen_t_Q2_%i",k);
      hgen_t_Q2[k]=new TH2F(hstname,hstname,200,0,2,250,0,20);
      sprintf(hstname,"hgen_theta_phi_CM_%i",k);
//    hgen_theta_phi_CM[k]=new TH2F(hstname,hstname,360,-3.15,3.15,180,0,3.15);      
      hgen_theta_phi_CM[k]=new TH2F(hstname,hstname,360,-3.15,3.15,100,-1,1);      
}


TFile *file=new TFile(input_filename.c_str());
    if (file->IsZombie()) {
       cout << "Error opening file" << input_filename << endl;
       return 0;
//        exit(-1);
}
else cout << "open file " << input_filename << endl;
TTree *Ttr1 = (TTree*) file->Get("tr1");

// Double_t Q2,t_t,psf,psf_flux,crs_BH,flux_factor;
// Double_t em_px,em_py,em_pz,ep_px,ep_py,ep_pz,prot_px,prot_py,prot_pz;
// Ttr1->SetBranchAddress("Q2",&Q2);
// Ttr1->SetBranchAddress("t",&t_t);
// Ttr1->SetBranchAddress("psf",&psf);
// Ttr1->SetBranchAddress("psf_flux",&psf_flux);
// Ttr1->SetBranchAddress("crs_BH",&crs_BH);
// Ttr1->SetBranchAddress("flux_factor",&flux_factor);
// Ttr1->SetBranchAddress("em_px",&em_px);
// Ttr1->SetBranchAddress("em_py",&em_py);
// Ttr1->SetBranchAddress("em_pz",&em_pz);
// Ttr1->SetBranchAddress("ep_px",&ep_px);
// Ttr1->SetBranchAddress("ep_py",&ep_py);
// Ttr1->SetBranchAddress("ep_pz",&ep_pz);
// Ttr1->SetBranchAddress("prot_px",&prot_px);
// Ttr1->SetBranchAddress("prot_py",&prot_py);
// Ttr1->SetBranchAddress("prot_pz",&prot_pz);

double Eg,Q2,t_t,psf,flux_factor,crs_BH,crs_INT;
TLorentzVector *L_em=0,*L_ep=0,*L_prot=0;
Ttr1->SetBranchAddress("Eg",&Eg);
Ttr1->SetBranchAddress("Q2",&Q2);
Ttr1->SetBranchAddress("t",&t_t);
Ttr1->SetBranchAddress("psf",&psf);
Ttr1->SetBranchAddress("flux_factor",&flux_factor);
Ttr1->SetBranchAddress("crs_BH",&crs_BH);
Ttr1->SetBranchAddress("crs_INT",&crs_INT);
Ttr1->SetBranchAddress("L_em",&L_em);
Ttr1->SetBranchAddress("L_ep",&L_ep);
Ttr1->SetBranchAddress("L_prot",&L_prot);


Int_t nevent = (Int_t)Ttr1->GetEntries();
Int_t nselected = 0;
cout << "nevent " << nevent << endl;

  Double_t cov= 1e-12 * 1e-24; //pb to cm2 conversion
  Double_t br = 1.;
  Double_t eff = 0.85;
  Double_t lumi, time;
if (detector=="SoLID"){
  lumi = 1.2e37;  // 1.2e37/cm2/s is from 3nA on 15cm long LH2 target
  time = 50*3600*24;  //50 days in seconds
}
else if (detector=="CLAS12"){
  lumi = 1e35;  // 1e35/cm2/s for CLAS12
  time = 130*3600*24;  //130 days in seconds
}
  Double_t overall_NOpsf = cov *lumi * br * eff * time /nevent;
  
  Int_t counter;
  Int_t countd = 1;
  cout << "|___________________|" << endl;
  cout << "|" << flush;

  int yescounter=0,nocounter=0;
  
      double weight[n]={0,0,0,0,0};
for (Int_t i=0;i<nevent;i++) {
 
    Ttr1->GetEntry(i);  

    cout << i << "\r"; 
    
    if (i<10)   cout << Eg << " " << Q2 << " " << t_t << " " << psf << " " << flux_factor  << crs_BH << " "  << crs_INT << " " << " " << L_em->M() << " " << L_ep->M() << " " << L_prot->M() << endl;
//     if (i<10)   cout << Q2 << " " << t_t << " " << psf << " " << psf_flux << " " << crs_BH << " " << flux_factor << " " << em_px<< " " <<em_py<< " " <<em_pz<< " " <<ep_px<< " " <<ep_py<< " " <<ep_pz<< " " <<prot_px<< " " <<prot_py<< " " <<prot_pz <<  endl;
    
    
    Double_t id = i;
    Double_t neventd = nevent;
    if(floor(id/neventd*20) == countd){   
      cout << "=" << flush;
      countd++;
      if(countd == 20){
	cout << "|" << endl;
      }
    }    
    
//   if (crs_BH > 1000000) {
//     cout << "crs > 1e6" << endl;
//       cout << Q2 << " " << t_t << " " << psf << " " << psf_flux << " " << crs_BH << " " << flux_factor << " " << em_px<< " " <<em_py<< " " <<em_pz<< " " <<ep_px<< " " <<ep_py<< " " <<ep_pz<< " " <<prot_px<< " " <<prot_py<< " " <<prot_pz <<  endl;
//   }

      hcrs_BH_log->Fill(log10(crs_BH));
  
//       if (isinf(crs_BH) || isnan(crs_BH) ) {cout << i << " " << crs_BH << endl; continue;}
//       if (isinf(crs_BH) || TMath::IsNaN(crs_BH) ) {nocounter++; continue;}
//       if (TMath::IsNaN(crs_BH)) {nocounter++; continue;}      
//      cout << i << " " << crs_BH << endl;
      
      if(Q2 < 4 || Q2 > 9) continue;            
      
      yescounter++;      

// 	 TLorentzVector *pr=L_prot,*e=L_em,*p=L_ep;
      TLorentzVector pr=*L_prot,e=*L_em,p=*L_ep;
//       TLorentzVector e,p,pr;
//       e.SetXYZM(em_px,em_py,em_pz,0.000511);
//       p.SetXYZM(ep_px,ep_py,ep_pz,0.000511);      
//       pr.SetXYZM(prot_px,prot_py,prot_pz,0.938);
      
      double pr_mom=pr.P(),e_mom=e.P(),p_mom=p.P();
      double pr_theta=pr.Theta(),e_theta=e.Theta(),p_theta=p.Theta();
      double pr_phi=pr.Phi(),e_phi=e.Phi(),p_phi=p.Phi();    
      
     
      ///get acc and resolution before smear
      int res_index_electron=0, res_index_positron=0, res_index_proton=0;
      double acc_electron_largeangle=0,acc_electron_forwardangle=0,acc_positron_largeangle=0,acc_positron_forwardangle=0,acc_proton_largeangle=0,acc_proton_forwardangle=0;
      double acc_proton=0,acc_positron=0,acc_electron=0;      
      bool cut=false;
      
      if (detector=="SoLID"){ 
       acc_proton=hacceptance_positive->GetBinContent(hacceptance_positive->FindBin(pr_theta*DEG,pr_mom));  
       acc_positron=hacceptance_positive->GetBinContent(hacceptance_positive->FindBin(p_theta*DEG,p_mom));
       acc_electron=hacceptance_negative->GetBinContent(hacceptance_negative->FindBin(e_theta*DEG,e_mom));    
       
       acc_electron_largeangle=hacceptance_negative_largeangle->GetBinContent(hacceptance_negative_largeangle->FindBin(e_theta*DEG,e_mom));
       acc_positron_largeangle=hacceptance_positive_largeangle->GetBinContent(hacceptance_positive_largeangle->FindBin(p_theta*DEG,p_mom));
       acc_proton_largeangle=hacceptance_positive_largeangle->GetBinContent(hacceptance_positive_largeangle->FindBin(pr_theta*DEG,pr_mom));
       
       acc_electron_forwardangle=hacceptance_negative_forwardangle->GetBinContent(hacceptance_negative_forwardangle->FindBin(e_theta*DEG,e_mom));
       acc_positron_forwardangle=hacceptance_positive_forwardangle->GetBinContent(hacceptance_positive_forwardangle->FindBin(p_theta*DEG,p_mom));
       acc_proton_forwardangle=hacceptance_positive_forwardangle->GetBinContent(hacceptance_positive_forwardangle->FindBin(pr_theta*DEG,pr_mom));
       
//        if ( (acc_electron_forwardangle>0 && e_mom < 4.9) || (acc_positron_forwardangle>0 && p_mom < 4.9) ) {}
//        else {acc_positron = 0 ; acc_electron = 0;}//cut on at least one lepton forward can use CC	 
//        if ( acc_electron_largeangle > 0 && acc_positron_largeangle > 0) {acc_positron=0; acc_electron=0;} //cut away two lepton in largeangle only   
//        if ( (acc_positron_largeangle>0 && acc_electron_forwardangle>0 && e_mom > 5.) || (acc_electron_largeangle>0 && acc_positron_forwardangle>0 && p_mom > 5.) ||  (acc_positron_forwardangle >0 && p_mom>5. && acc_electron_forwardangle>0 && e_mom > 5.) ) {acc_positron=0; acc_electron=0;} //cut away forward lepton exceed CC limit       
       if (acc_proton_largeangle > 0 && pr_mom > 2.0)  acc_proton=0; //cut away 2.0 at largeangle  for proton
       if (acc_proton_forwardangle > 0 && pr_mom > 4.5)  acc_proton=0; //cut away 4.4 at forwardangle for proton
       
       if (acc_electron_forwardangle>0) res_index_electron=0;
       if (acc_electron_largeangle  >0) res_index_electron=1;
       if (acc_positron_forwardangle>0) res_index_positron=0;
       if (acc_positron_largeangle  >0) res_index_positron=1;
       if (acc_proton_forwardangle  >0) res_index_proton=0;
       if (acc_proton_largeangle    >0) res_index_proton=1;
      }   
      else if (detector=="CLAS12"){
	  //note CLAS12 phi (0,360), TCS sim phi (-180,180) 	
	  acc_proton=hacceptance_PThetaPhi_positive->GetBinContent(hacceptance_PThetaPhi_positive->FindBin(pr_phi*DEG+180,pr_theta*DEG,pr_mom));
	  acc_positron=hacceptance_PThetaPhi_positive->GetBinContent(hacceptance_PThetaPhi_positive->FindBin(p_phi*DEG+180,p_theta*DEG,p_mom));
	  acc_electron=hacceptance_PThetaPhi_negative->GetBinContent(hacceptance_PThetaPhi_negative->FindBin(e_phi*DEG+180,e_theta*DEG,e_mom));	  
	  
	  if (p_theta*DEG>36) acc_positron = 0; //lepton has no detection at central detector
	  if (e_theta*DEG>36) acc_electron = 0; //lepton has no detection at central detector
	  if (e_mom > 4.9 && p_mom > 4.9) {acc_positron=0;acc_electron=0;} // at least one in CC limit
	  if (pr_theta*DEG>36 && pr_mom > 1) acc_proton=0;
	  if (pr_theta*DEG<36 && pr_mom > 4) acc_proton=0;	  
	  //cut away unsure low acceptance	  	  
	  if (acc_proton < 0.9) acc_proton=0;
	  if (acc_positron < 0.9) acc_positron=0;
	  if (acc_electron < 0.9) acc_electron=0; 
	  
	  if (e_theta*DEG<=36) res_index_electron=0;
	  if (e_theta*DEG >36) res_index_electron=1;
	  if (p_theta*DEG<=36) res_index_positron=0;
	  if (p_theta*DEG >36) res_index_positron=1;
	  if (pr_theta*DEG<=36) res_index_proton=0;
	  if (pr_theta*DEG >36) res_index_proton=1;
      }
      else if (detector=="EIC"){
      }
      else {return 0;}
   
      //for trigger study
//       if (detector=="SoLID"){ 
//        acc_proton=hacceptance_positive->GetBinContent(hacceptance_positive->FindBin(pr_theta*DEG,pr_mom));  
//        acc_positron=hacceptance_positive->GetBinContent(hacceptance_positive->FindBin(p_theta*DEG,p_mom));
//        acc_electron=hacceptance_negative->GetBinContent(hacceptance_negative->FindBin(e_theta*DEG,e_mom));    
//        
//        acc_electron_largeangle=hacceptance_negative_largeangle->GetBinContent(hacceptance_negative_largeangle->FindBin(e_theta*DEG,e_mom));
//        acc_positron_largeangle=hacceptance_positive_largeangle->GetBinContent(hacceptance_positive_largeangle->FindBin(p_theta*DEG,p_mom));
//        acc_proton_largeangle=hacceptance_positive_largeangle->GetBinContent(hacceptance_positive_largeangle->FindBin(pr_theta*DEG,pr_mom));
//        
//        acc_electron_forwardangle=hacceptance_negative_forwardangle->GetBinContent(hacceptance_negative_forwardangle->FindBin(e_theta*DEG,e_mom));
//        acc_positron_forwardangle=hacceptance_positive_forwardangle->GetBinContent(hacceptance_positive_forwardangle->FindBin(p_theta*DEG,p_mom));
//        acc_proton_forwardangle=hacceptance_positive_forwardangle->GetBinContent(hacceptance_positive_forwardangle->FindBin(pr_theta*DEG,pr_mom));
//        
//        if (e_mom < 2 || p_mom < 2)  cut=true;
//        
//        if (acc_electron_forwardangle>0) res_index_electron=0;
//        if (acc_electron_largeangle  >0) res_index_electron=1;
//        if (acc_positron_forwardangle>0) res_index_positron=0;
//        if (acc_positron_largeangle  >0) res_index_positron=1;
//        if (acc_proton_forwardangle  >0) res_index_proton=0;
//        if (acc_proton_largeangle    >0) res_index_proton=1;
//       } 
      
      double acc=acc_proton*acc_positron*acc_electron;         
//       double acc=acc_positron*acc_electron;               
     
      ///smear by detector resolution 
      if(Is_smear){
      double temp_p,temp_theta,temp_phi;
      temp_p=gRandom->Gaus(pr.P(),pr.P()*resolution_p[res_index_proton]);
      temp_theta=gRandom->Gaus(pr.Theta(),resolution_theta[res_index_proton]);
      temp_phi=gRandom->Gaus(pr.Phi(),resolution_theta[res_index_proton]);
      pr.SetXYZM(temp_p*sin(temp_theta)*cos(temp_phi),temp_p*sin(temp_theta)*sin(temp_phi),temp_p*cos(temp_theta),0.938);
      temp_p=gRandom->Gaus(p.P(),p.P()*resolution_p[res_index_positron]);
      temp_theta=gRandom->Gaus(p.Theta(),resolution_theta[res_index_positron]);
      temp_phi=gRandom->Gaus(p.Phi(),resolution_theta[res_index_positron]);
      p.SetXYZM(temp_p*sin(temp_theta)*cos(temp_phi),temp_p*sin(temp_theta)*sin(temp_phi),temp_p*cos(temp_theta),0.000511);
      temp_p=gRandom->Gaus(e.P(),e.P()*resolution_p[res_index_electron]);
      temp_theta=gRandom->Gaus(e.Theta(),resolution_theta[res_index_electron]);
      temp_phi=gRandom->Gaus(e.Phi(),resolution_theta[res_index_electron]);
      e.SetXYZM(temp_p*sin(temp_theta)*cos(temp_phi),temp_p*sin(temp_theta)*sin(temp_phi),temp_p*cos(temp_theta),0.000511);    
      }
  
      double M=0.938;      
//       TLorentzVector target(0.,0.,0.,M);                  
      TLorentzVector ph=pr+e+p-target;
      TLorentzVector ep=e+p;
      TLorentzVector beam(0.,0.,11.,11.);
      TLorentzVector Miss=beam-ph;
      
//       cout << ph.P() << " " << ph.E() << " " << ph.Theta() << " " << ph.Phi()  << endl;
//       ph.SetXYZM(0,0,ph.P(),0);
   
      pr_mom=pr.P();e_mom=e.P();p_mom=p.P();
      pr_theta=pr.Theta();e_theta=e.Theta();p_theta=p.Theta();
      pr_phi=pr.Phi();e_phi=e.Phi();p_phi=p.Phi();
      
      double ph_mom=ph.P(),ph_theta=ph.Theta(),ph_phi=ph.Phi();
      
      double InvM_ep=ep.M();
      double t=-(pr-target).M2();
      double s=(ph+target).M2();
      double W=(ph+target).M();
      double tau=InvM_ep*InvM_ep/(s-M*M);
//       double eta=tau/(2-tau);   
      double eta=Q2/(2*s-2*M*M-Q2);   //       tau/(2-tau)=Q2/(2*s-2*M*M-Q2)
//       double eta=Q2/(2*(s-M*M)-Q2-t);  //this is detailed defination  
      double xi=Q2/(2*s-Q2);       //eta in TCS is similar to xi in DVCS 
      double t_min=4*eta*eta*M*M/(1-eta*eta);   //???????????????      
     
//       cout << " t " << t_t << " " << -t << endl;
//       cout << " M " << Q2 << " " << InvM_ep*InvM_ep << endl;      
//       cout << InvM_ep << " " << t << " " << s << " " << W << endl; 

      TLorentzVector target_CM=target,pr_CM=pr,e_CM=e,p_CM=p; 
      target_CM.Boost(-ep.BoostVector()); 
      pr_CM.Boost(-ep.BoostVector());
      e_CM.Boost(-ep.BoostVector());
      p_CM.Boost(-ep.BoostVector());
      TVector3 e3=-pr_CM.Vect().Unit();
      TVector3 e2=(e3.Cross(target_CM.Vect())).Unit();
      TVector3 e1=e2.Cross(e3);                          
      double z=e_CM.Vect().Dot(e3),y=e_CM.Vect().Dot(e2),x=e_CM.Vect().Dot(e1);
      double theta_CM=acos(z/sqrt(x*x+y*y+z*z));
      double phi_CM=acos(x/sqrt(x*x+y*y));
      if (y<0) phi_CM=-phi_CM;
//       cout << " phi_CM " << phi_cm  << " " << phi_CM << endl;
//       cout << " th_cm " << th_cm << " " << theta_CM << endl;   

      ///================verify inverse derive the kinematics =========
//       double m=0.000511;
//       double Pi=3.1415926;
//       double       Eg=(s/M-M)/2;      
//       double       prot_E=t/2/M+M;
//       double       prot_P=sqrt(prot_E*prot_E-M*M);
//       double       prot_Phi=pr_phi;      
//       double       tp_E=Eg+M-prot_E;
//       double       tp_P=sqrt(tp_E*tp_E-Q2);
//       double       tp_Phi=prot_Phi-Pi*prot_Phi/fabs(prot_Phi);      
//       double       prot_Theta=acos((Eg*Eg-tp_P*tp_P+prot_P*prot_P)/(2*Eg*prot_P));            
//       double       tp_Theta=acos((Eg*Eg+tp_P*tp_P-prot_P*prot_P)/(2*Eg*tp_P));      
//       double       ele_CM_E=sqrt(Q2)/2;
//       double       ele_CM_P=sqrt(ele_CM_E*ele_CM_E-m*m);
//       double       ele_CM_Theta=theta_CM;
//       double       ele_CM_Phi=phi_CM;      
//       double       pos_CM_E=ele_CM_E;
//       double       pos_CM_P=ele_CM_P;
//       double       pos_CM_Theta=Pi-ele_CM_Theta;
//       double       pos_CM_Phi=ele_CM_Phi-Pi*ele_CM_Phi/fabs(ele_CM_Phi); 
//       TVector3  ele_CM_v=ele_CM_P*sin(ele_CM_Theta)*cos(ele_CM_Phi)*e1+ele_CM_P*sin(ele_CM_Theta)*sin(ele_CM_Phi)*e2+ele_CM_P*cos(ele_CM_Theta)*e3;
//       TVector3  pos_CM_v=pos_CM_P*sin(pos_CM_Theta)*cos(pos_CM_Phi)*e1+pos_CM_P*sin(pos_CM_Theta)*sin(pos_CM_Phi)*e2+pos_CM_P*cos(pos_CM_Theta)*e3;
//       TLorentzVector ele_CM(ele_CM_v,ele_CM_E); 
//       TLorentzVector pos_CM(pos_CM_v,pos_CM_E);
//       TLorentzVector tp(tp_P*sin(tp_Theta)*cos(tp_Phi),tp_P*sin(tp_Theta)*sin(tp_Phi),tp_P*cos(tp_Theta),tp_E);      
//       TLorentzVector ele=ele_CM,pos=pos_CM;
//       ele.Boost(tp.BoostVector());
//       pos.Boost(tp.BoostVector());           
//       double       ele_E=ele.E();
//       double       ele_P=ele.P();
//       double       ele_Theta=ele.Theta();
//       double       ele_Phi=ele.Phi();
//       double       pos_E=pos.E();
//       double       pos_P=pos.P();
//       double       pos_Theta=pos.Theta();
//       double       pos_Phi=pos.Phi();

//       cout << pr_mom << " " <<  pr_theta << " " <<  pr_phi << "\t " << ep.P() << " " <<  ep.Theta() << " " <<  ep.Phi() << "\t " << e_mom<< " " << e_theta<< " " << e_phi<< "\t " << p_mom<< " " << p_theta<< " " << p_phi << endl;    
//       cout << prot_P << " " <<  prot_Theta << " " <<  prot_Phi << "\t " << tp_P << " " <<  tp_Theta << " " <<  tp_Phi  << "\t " <<  ele_P << " " <<  ele_Theta << " " <<  ele_Phi  << "\t " << pos_P << " " <<  pos_Theta << " " <<  pos_Phi << endl;
      ///====================================================================
      
      
      
      hEg->Fill(ph_mom);
      hflux_factor->Fill(flux_factor);
      hEgflux_factor->Fill(ph_mom,flux_factor);
      hEg_flux_factor->Fill(ph_mom,flux_factor);      

      // cut to quasi real
      double acc_cut=1;
//       if (fabs(ph.M2())>0.02) acc_cut=0;
      if (fabs(ph.M2())>0.05) acc_cut=0;      
//       if(fabs(Miss.M2())>0.05 || fabs(Miss.Px()/Miss.P())>0.2 || fabs(Miss.Py()/Miss.P())>0.2) acc_cut=0;      

// 	cout << acc << " " << acc_proton << " " << acc_positron << " " << acc_electron << endl;
// 	cout << pr_theta*DEG << " " << pr_phi*DEG << " " << p_theta*DEG << " " << p_phi*DEG << " " << e_theta*DEG << " " << e_phi*DEG << endl;
  
      weight[0]=1;
//       weight[1]=crs_BH*psf*psf_flux*flux_factor;
      weight[1]=crs_BH*psf*flux_factor;      
      weight[2]=weight[1]*acc;
      weight[3]=weight[2]*acc_cut;
      weight[4]=weight[3]*overall_NOpsf;      
//       weight[3]=weight[2]*overall_NOpsf;
//       weight[4]=weight[3]*acc_cut;
//       double weight[n]={1,crs_BH*psf*psf_flux*flux_factor,crs_BH*psf*psf_flux*flux_factor*acc,crs_BH*psf*psf_flux*flux_factor*acc*overall_NOpsf};
      for(int k=0;k<n;k++){
	
	heta_xi[k]->Fill(xi,eta,weight[k]);
	
	ht_Q2_xi[k]->Fill(xi,Q2,t,weight[k]);
	hQ2_xi[k]->Fill(xi,Q2,weight[k]);
  
	hphoton_m2[k]->Fill(fabs(ph.M2()),weight[k]);
	
	hphoton_theta_mom[k]->Fill(ph_theta*DEG,ph_mom,weight[k]);
	hproton_theta_mom[k]->Fill(pr_theta*DEG,pr_mom,weight[k]);
	helectron_theta_mom[k]->Fill(e_theta*DEG,e_mom,weight[k]);
	hpositron_theta_mom[k]->Fill(p_theta*DEG,p_mom,weight[k]);
	
	hep_mom[k]->Fill(p_mom,e_mom,weight[k]);
	hep_theta[k]->Fill(p_theta*DEG,e_theta*DEG,weight[k]);
	ht[k]->Fill(t,weight[k]);
	hs[k]->Fill(s,weight[k]);
	htau[k]->Fill(tau,weight[k]);
	heta[k]->Fill(eta,weight[k]);
	hep_InvM[k]->Fill(InvM_ep,weight[k]);
	hQ2_t[k]->Fill(t,Q2,weight[k]);
	htheta_phi_CM[k]->Fill(phi_CM,theta_CM,weight[k]);
	hphi_CM[k]->Fill(phi_CM,weight[k]);
	
        hMissP_MM2[k]->Fill(Miss.M2(),Miss.P(),weight[k]);
        hQ2_MM2[k]->Fill(Miss.M2(),-ph.M2(),weight[k]);	
        hMissPxPy[k]->Fill(Miss.Px()/Miss.P(),Miss.Py()/Miss.P(),weight[k]);
	hMissMM2[k]->Fill(Miss.M2(),weight[k]);
        hMissPt[k]->Fill(Miss.Pt()/Miss.P(),weight[k]);
	
	if (17.5 < s && s < 19.5){
	  if (0.1 < t && t < 0.9){
	    htheta_phi_CM_bin[int((t-0.1)/0.2)][k]->Fill(phi_CM,theta_CM,weight[k]);
  
// 	    cout << InvM_ep << " " << s << endl;
// 	    if (3.1415926/4.<theta_CM && theta_CM <3.*3.1415926/4. && fabs(InvM_ep*InvM_ep-5.)<0.1 && fabs(sqrt(s)-sqrt(20))<0.1){	      
// 		hphi_CM_bin[int((t-0.1)/0.2)][k]->Fill(phi_CM,weight[k]);
// 	    }
	    
	    hphi_CM_bin[int((t-0.1)/0.2)][k]->Fill(phi_CM,weight[k]);
	  }
	}
	
      }
      
	if (17.5 < s && s < 19.5){
	  if (6 < Q2 && Q2 < 7){
// 	    hcrs_t_stat->Fill(t,crs_BH*psf*psf_flux*flux_factor/nevent*acc*acc_cut);
// 	    hcrs_t->Fill(t,crs_BH*psf*psf_flux*flux_factor/nevent/0.2);
	    hcrs_t_stat->Fill(t,crs_BH*psf*flux_factor/nevent*acc*acc_cut);
	    hcrs_t->Fill(t,crs_BH*psf*flux_factor/nevent/0.2);	    
	  }
	}
	
	if ((24 < s && s < 26) && (4.5 < Q2 && Q2 < 5.5) && (3.1416/4 < theta_CM && theta_CM < 3*3.1416/4)){
// 	if ((0 < s && s < 30) && (6.4 < Q2 && Q2 < 6.6) && (3.1416/4. < theta_CM && theta_CM < 3.*3.1416/4.)){  
	    hcrs_t_paper->Fill(t,crs_BH*psf/nevent/0.01/0.2);
	    hcrs_t_paper_nopsf->Fill(t,crs_BH/nevent/0.01/0.2);
	}
  
	
	hs_final->Fill(s,weight[index]);	
	ht_final->Fill(t,weight[index]);
	heta_final->Fill(eta,weight[index]);
	hQ2_final->Fill(Q2,weight[index]);
	hs_t_final->Fill(s,t,weight[index]);
	hs_Q2_final->Fill(s,Q2,weight[index]);
	heta_t_final->Fill(eta,t,weight[index]);
	heta_Q2_final->Fill(eta,Q2,weight[index]);
	ht_Q2_final->Fill(t,Q2,weight[index]);
	heta_t_Q2_final->Fill(eta,t,Q2,weight[index]);
	hs_t_Q2_final->Fill(s,t,Q2,weight[index]);	

	hs_final_bin->Fill(s,weight[index]);	
	ht_final_bin->Fill(t,weight[index]);
	heta_final_bin->Fill(eta,weight[index]);
	hQ2_final_bin->Fill(Q2,weight[index]);
	hs_t_final_bin->Fill(s,t,weight[index]);
	hs_Q2_final_bin->Fill(s,Q2,weight[index]);
	heta_t_final_bin->Fill(eta,t,weight[index]);
	heta_Q2_final_bin->Fill(eta,Q2,weight[index]);
	ht_Q2_final_bin->Fill(t,Q2,weight[index]);
	heta_t_Q2_final_bin->Fill(eta,t,Q2,weight[index]);
	hs_t_Q2_final_bin->Fill(s,t,Q2,weight[index]);		
        
      if ((t-t_min) > 0.2) continue;
      
      int sbin_even=int((s-10)/((25-10)/Nsbin));
      int etabin_even=int((eta-0.1)/((0.4-0.1)/Netabin));
      int Q2bin_even=int((Q2-4.)/((9.-4.)/NQ2bin));
      int tbin_even=int((t-0.)/((4.-0.)/Ntbin));
      
      if ( (0 <= etabin_even && etabin_even < Netabin) && (0 <= Q2bin_even && Q2bin_even < NQ2bin) && (0<= tbin_even && tbin_even < Ntbin) ){  
	ht_Q2_etabin_even[etabin_even]->Fill(t,Q2,weight[index]);      
      }
      
      if ( (0 <= sbin_even && sbin_even < Nsbin) && (0 <= Q2bin_even && Q2bin_even < NQ2bin) && (0<= tbin_even && tbin_even < Ntbin) ){  
	ht_Q2_sbin_even[sbin_even]->Fill(t,Q2,weight[index]);      
      }     
	
	int sbin=-1,etabin=-1,Q2bin=-1,tbin=-1;
	for (int k=0;k<Nsbin;k++) if (*(sbin_edge+k) <= s && s < *(sbin_edge+k+1)) sbin=k;	
	for (int k=0;k<Netabin;k++) if (*(etabin_edge+k) <= eta && eta < *(etabin_edge+k+1)) etabin=k;
	for (int k=0;k<NQ2bin;k++)  if (*(Q2bin_edge+k) <= Q2 && Q2 < *(Q2bin_edge+k+1)) Q2bin=k;
	for (int k=0;k<Ntbin;k++) if (*(tbin_edge+k) <= t && t < *(tbin_edge+k+1)) tbin=k;
	
	if (0 <= etabin && etabin < Netabin){
	    ht_Q2_etabin[etabin]->Fill(t,Q2,weight[index]);
	    htheta_phi_CM_etabin[etabin]->Fill(phi_CM,theta_CM,weight[index]);
	    int Q2bin_etabin=-1;
	    for (int k=0;k<NQ2bin_etabin;k++)  if (*(Q2bin_etabin_edge+etabin*(NQ2bin_etabin+1)+k) <= Q2 && Q2 < *(Q2bin_etabin_edge+etabin*(NQ2bin_etabin+1)+k+1)) Q2bin_etabin=k;
	    if ( (0 <= Q2bin_etabin && Q2bin_etabin < NQ2bin_etabin) ){
	      countbin_etabin[etabin][Q2bin_etabin][0] += weight[index];
	    }  
	}
	
	if (0 <= Q2bin && Q2bin < NQ2bin){
	    ht_eta_Q2bin[Q2bin]->Fill(t,eta,weight[index]);
	    htheta_phi_CM_Q2bin[Q2bin]->Fill(phi_CM,theta_CM,weight[index]);
	    int etabin_Q2bin=-1;
	    for (int k=0;k<Netabin_Q2bin;k++)  if (*(etabin_Q2bin_edge+Q2bin*(Netabin_Q2bin+1)+k) <= eta && eta < *(etabin_Q2bin_edge+Q2bin*(Netabin_Q2bin+1)+k+1)) etabin_Q2bin=k;
	    if ( (0 <= etabin_Q2bin && etabin_Q2bin < Netabin_Q2bin) ){
	      countbin_Q2bin[Q2bin][etabin_Q2bin][0] += weight[index];
	    }  
	}	
	
//       if ( (0 <= etabin && etabin < Netabin) && (0 <= Q2bin && Q2bin < NQ2bin) && (0<= tbin && tbin < Ntbin) ){  
// 	cout << "in range " << etabin << " " << tbin << " " << Q2bin << endl;
// 	hThetaPhiCM_etabin_Q2bin_tbin[etabin][Q2bin][tbin]->Fill(phi_CM,theta_CM,weight[index]);
// 	vector<TH2F>::iterator iter=hThetaPhiCM_etabin_Q2bin_tbin->begin()+(etabin*Netabin+tbin*Ntbin+Q2bin*NQ2bin);
// 	iter->Fill(phi_CM,theta_CM,weight[n-1]);

// 	countbin_etabin[etabin][Q2bin][tbin] = countbin_etabin[etabin][Q2bin][tbin] + weight[index];
//       }
//       else {
// 	cout << etabin << " " << Q2bin  << " " << tbin << endl;
// 	cout << "out of range " << eta << " " << Q2 << " " << t << endl;	
// 	nocounter++;
//       } 
      
//       if ( (0 <= sbin && sbin < Nsbin) && (0 <= Q2bin && Q2bin < NQ2bin) && (0<= tbin && tbin < Ntbin) ){  
      if ( (0 <= sbin && sbin < Nsbin) ){    
	ht_Q2_sbin[sbin]->Fill(t,Q2,weight[index]);
      }
      else {
// 	cout << sbin << " " << Q2bin  << " " << tbin << endl;
// 	cout << "out of range " << s << " " << Q2 << " " << t << endl;	
	nocounter++;
      }       

            	
//       if (6 < ph_mom && ph_mom < 11) {
// 	hgen_t_Q2[int((ph_mom-6)/0.25)]->Fill(t,Q2,psf);
// 	hgen_t_Q2[int((ph_mom-6)/0.25)]->Fill(t,InvM_ep);      
      
// 	if (0.5<t && t <1.6 && 4 < Q2 && Q2< 9) hgen_theta_phi_CM[int((ph_mom-0.)/0.25)]->Fill(phi_CM,cos(theta_CM));
//       }

}

      cout << "yescounter " << yescounter << " nocounter " << nocounter << endl;      

// for(int etabin=0;etabin<Netabin;etabin++){
// cout << "etabin " << etabin << " " <<  Q2_min_etabin[etabin] << " " << Q2_max_etabin[etabin] << " " << t_min_etabin[etabin] << " " << t_max_etabin[etabin] << endl;
// for(int Q2bin=0;Q2bin<5;Q2bin++){
// for(int tbin=0;tbin<5;tbin++){  
// //   double a=hThetaPhiCM_etabin_Q2bin_tbin[etabin][Q2bin][tbin]->Integral();
// //   cout << etabin << " " << Q2bin  << " " << tbin << " " << a << endl;
//   
// //   cout << etabin << " " << tbin << " " << Q2bin << endl;
// //   Double_t a=countbin_etabin[etabin][Q2bin][tbin];
// //   if (isinf(a) || TMath::IsNaN(a)) cout << etabin << " " << tbin << " " << Q2bin << " " << " bad" <<  endl;
// 
//   cout << etabin << " " << Q2bin << " " <<  tbin << " "  << countbin_etabin[etabin][Q2bin][tbin] <<  endl;
// 
//   //     vector<TH2F>::iterator iter=hThetaPhiCM_etabin_Q2bin_tbin->begin()+(etabin*Netabin+tbin*Ntbin+Q2bin*NQ2bin);
// //     cout << etabin << " " << tbin << " " << Q2bin << " " << iter->Integral() <<  endl;    
//       
// }}}        

TCanvas *c_theta_mom = new TCanvas("theta_mom","theta_mom",1500,900);
c_theta_mom->Divide(4,n);
for(int k=0;k<n;k++){
//   double max,min;
//   max=log10(hproton_theta_mom[k]->GetBinContent(hproton_theta_mom[k]->GetMaximumBin()))+1;
//   min=log10(hproton_theta_mom[k]->GetBinContent(hproton_theta_mom[k]->GetMinimumBin()))-1;
  
  c_theta_mom->cd(4*k+1);
  gPad->SetLogz(1);    
  hphoton_theta_mom[k]->Draw("colz");
  c_theta_mom->cd(4*k+2);
  gPad->SetLogz(1);  
//  hproton_theta_mom[k]->SetMinimum(min);
 // hproton_theta_mom[k]->SetMaximum(max);     
//   hproton_theta_mom[k]->SetAxisRange(0,60,"X");
  hproton_theta_mom[k]->Draw("colz");
  c_theta_mom->cd(4*k+3);
  gPad->SetLogz(1);
//   helectron_theta_mom[k]->SetMinimum(min);
//   helectron_theta_mom[k]->SetMaximum(max);
//   helectron_theta_mom[k]->SetAxisRange(0,40,"X");  
  helectron_theta_mom[k]->Draw("colz");
  c_theta_mom->cd(4*k+4);
  gPad->SetLogz(1);
// hpositron_theta_mom[k]->SetMinimum(min);
 // hpositron_theta_mom[k]->SetMaximum(max);  
//   hpositron_theta_mom[k]->SetAxisRange(0,40,"X");  
  hpositron_theta_mom[k]->Draw("colz");
//   c_theta_mom->cd(5*k+5);
//   gPad->SetLogz(1);  
//   helectron_positron_theta_mom_ratio[k]->Divide(helectron_theta_mom[k],hpositron_theta_mom[k]);  
// //   max=log10(hproton_theta_mom[k]->GetBinContent(hproton_theta_mom[k]->GetMaximumBin()))+1;
// //   min=log10(hproton_theta_mom[k]->GetBinContent(hproton_theta_mom[k]->GetMinimumBin()))-1;
// //   helectron_positron_theta_mom_ratio[k]->SetMinimum(min);
// //   helectron_positron_theta_mom_ratio[k]->SetMaximum(max);
// //   helectron_positron_theta_mom_ratio[k]->SetAxisRange(0,40,"X");
//   helectron_positron_theta_mom_ratio[k]->Draw("colz");  
}

TCanvas *c_theta_mom_final = new TCanvas("theta_mom_final","theta_mom_final",1200,1000);
c_theta_mom_final->Divide(2,2);
for(int k=index;k<index+1;k++){
  c_theta_mom_final->cd(1);
  gPad->SetLogz(1); 
  hphoton_theta_mom[k]->SetMinimum(0.01);
  hphoton_theta_mom[k]->SetAxisRange(0,thetamax,"X");  
  hphoton_theta_mom[k]->Draw("colz");
  c_theta_mom_final->cd(2);
  gPad->SetLogz(1);
  hproton_theta_mom[k]->SetMinimum(0.01);
  hproton_theta_mom[k]->SetAxisRange(0,thetamax,"X");  
  hproton_theta_mom[k]->Draw("colz");
  c_theta_mom_final->cd(3);
  gPad->SetLogz(1);
  helectron_theta_mom[k]->SetMinimum(0.01);
  helectron_theta_mom[k]->SetAxisRange(0,thetamax,"X");    
  helectron_theta_mom[k]->Draw("colz");
  c_theta_mom_final->cd(4);
  gPad->SetLogz(1);
  hpositron_theta_mom[k]->SetMinimum(0.01);
  hpositron_theta_mom[k]->SetAxisRange(0,thetamax,"X");    
  hpositron_theta_mom[k]->Draw("colz");
}

TCanvas *c_ep_final = new TCanvas("ep_final","ep_final",1200,600);
c_ep_final->Divide(2,1);
c_ep_final->cd(1);
gPad->SetLogz();
// hep_mom[index]->SetMinimum(0.01);
hep_mom[index]->Draw("colz");
c_ep_final->cd(2);
// hep_theta[index]->SetMinimum(0.01);
gPad->SetLogz();
hep_theta[index]->SetAxisRange(0,40,"X");
hep_theta[index]->SetAxisRange(0,40,"Y");
hep_theta[index]->Draw("colz");

int nc=11;
TCanvas *c_other = new TCanvas("other","other",1800,900);
c_other->Divide(nc,n);
for(int k=0;k<n;k++){
  c_other->cd(nc*k+1);
  gPad->SetLogz(1);
  hep_mom[k]->SetMinimum(0.01);
  hep_mom[k]->Draw("colz");
  c_other->cd(nc*k+2);
  gPad->SetLogz(1);
  hep_theta[k]->SetMinimum(0.01);
  hep_theta[k]->Draw("colz");
  c_other->cd(nc*k+3);
  gPad->SetLogy(1);  
  ht[k]->Draw();
  c_other->cd(nc*k+4);
  gPad->SetLogy(1);  
  hs[k]->Draw();
  c_other->cd(nc*k+5);
  gPad->SetLogy(1);  
  hep_InvM[k]->Draw();
  c_other->cd(nc*k+6);
  gPad->SetLogz(1);
  hQ2_t[k]->Draw("colz");
  c_other->cd(nc*k+7);
  gPad->SetLogy(1);  
  htau[k]->Draw();
  c_other->cd(nc*k+8);
  gPad->SetLogy(1);  
  heta[k]->Draw();
  c_other->cd(nc*k+9);
  gPad->SetLogz(1);
  htheta_phi_CM[k]->SetMinimum(0.01);
  htheta_phi_CM[k]->Draw("colz");
  c_other->cd(nc*k+10);
  hphi_CM[k]->Draw();
  c_other->cd(nc*k+11);
  gPad->SetLogx(1);
//   gPad->SetLogy(1); 
  hphoton_m2[k]->Draw();  
  
}

TCanvas *c_theta_phi_CM_bin = new TCanvas("theta_phi_CM_bin","theta_phi_CM_bin",1600,900);
c_theta_phi_CM_bin->Divide(4,n);
for(int k=0;k<n;k++){
  for(int l=0;l<4;l++){
   c_theta_phi_CM_bin->cd(4*k+l+1);
   gPad->SetLogz(1);
   htheta_phi_CM_bin[l][k]->SetMinimum(0.01);   
   htheta_phi_CM_bin[l][k]->Draw("colz");
  }
}

TCanvas *c_theta_phi_CM_bin_final = new TCanvas("theta_phi_CM_bin_final","theta_phi_CM_bin_final",1000,800);
c_theta_phi_CM_bin_final->Divide(2,2);
for(int l=0;l<4;l++){
   c_theta_phi_CM_bin_final->cd(l+1);
   gPad->SetLogz(1);
   htheta_phi_CM_bin[l][3]->SetMinimum(0.01);   
   htheta_phi_CM_bin[l][3]->Draw("colz");  
   cout << htheta_phi_CM_bin[l][3]->Integral() << endl;  
}

TCanvas *c_phi_CM_bin = new TCanvas("phi_CM_bin","phi_CM_bin",1600,900);
c_phi_CM_bin->Divide(4,n);
for(int k=0;k<n;k++){
  for(int l=0;l<4;l++){
   c_phi_CM_bin->cd(4*k+l+1);
   hphi_CM_bin[l][k]->Draw();
  }
}

TCanvas *c_gen_t_Q2 = new TCanvas("gen_t_Q2","gen_t_Q2",1600,900);
c_gen_t_Q2->Divide(5,4);
for(int k=0;k<20;k++){
   c_gen_t_Q2->cd(k+1);
   hgen_t_Q2[k]->Draw("colz");
}

TCanvas *c_gen_theta_phi_CM = new TCanvas("gen_theta_phi_CM","gen_theta_phi_CM",1600,900);
c_gen_theta_phi_CM->Divide(5,4);
for(int k=0;k<20;k++){
   c_gen_theta_phi_CM->cd(k+1);
   hgen_theta_phi_CM[k]->Draw("colz");
}

TCanvas *c_crs_BH_log = new TCanvas("crs_BH_log","crs_BH_log",600,600);
gPad->SetLogy(1);
hcrs_BH_log->Draw();


cout << "crs_t error" << endl;
for(int k=1;k<6;k++){
// double err_rel=1./sqrt(hcrs_t->GetBinContent(k)/50.);
// double crs = hcrs_t->GetBinContent(k)/(50*3600*24*1e-36*1e37*0.85)/0.2;
double crs = hcrs_t->GetBinContent(k);
double err_rel=1./sqrt(hcrs_t_stat->GetBinContent(k)*1e-36*1e35*100*3600*24*0.85);
double err = crs*err_rel;
cout << hcrs_t->GetBinContent(k) << " " << err_rel << " " << err << " " << crs << endl;
hcrs_t->SetBinContent(k,crs);
hcrs_t->SetBinError(k,err);
}
TCanvas *c_crs_t = new TCanvas("crs_t","crs_t",600,600);
// gPad->SetLogy(0);
hcrs_t->SetMarkerStyle(24);
hcrs_t->SetMarkerColor(kBlack);
hcrs_t->SetTitle(";t(GeV2);d#sigma/dt (pb/GeV2)");
hcrs_t->Draw("E");

TCanvas *c_crs_t_paper = new TCanvas("crs_t_paper","crs_t_paper",600,600);
// gPad->SetLogy();
hcrs_t_paper->SetMarkerStyle(24);
hcrs_t_paper->SetMarkerColor(kBlack);
hcrs_t_paper->SetTitle(";t(GeV2);d#sigma/dt/dQ'2 (pb/GeV4)");
hcrs_t_paper->Draw("E");
hcrs_t_paper_nopsf->SetMarkerStyle(24);
hcrs_t_paper_nopsf->SetMarkerColor(kRed);
hcrs_t_paper_nopsf->SetTitle(";t(GeV2);d#sigma/dt/dQ'2 (pb/GeV4)");
// hcrs_t_paper_nopsf->Draw("E same");


TCanvas *c_Miss_2D = new TCanvas("Miss_2D","Miss_2D",1300,600);
c_Miss_2D->Divide(3,2);
c_Miss_2D->cd(1);
gPad->SetLogz();
hQ2_MM2[index-1]->Draw("colz");
c_Miss_2D->cd(2);
gPad->SetLogz();
hMissP_MM2[index-1]->Draw("colz");
c_Miss_2D->cd(3);
gPad->SetLogz();
hMissPxPy[index-1]->Draw("colz");
c_Miss_2D->cd(4);
gPad->SetLogz();
hQ2_MM2[index]->Draw("colz");
c_Miss_2D->cd(5);
gPad->SetLogz();
hMissP_MM2[index]->Draw("colz");
c_Miss_2D->cd(6);
gPad->SetLogz();
hMissPxPy[index]->Draw("colz");

TCanvas *c_Miss_1D = new TCanvas("Miss_1D","Miss_1D",1800,600);
c_Miss_1D->Divide(3,1);
c_Miss_1D->cd(1);
gPad->SetLogx();
hphoton_m2[index-1]->Draw();
c_Miss_1D->cd(2);
gPad->SetLogy();
hMissMM2[index-1]->Draw();
c_Miss_1D->cd(3);
gPad->SetLogy();
hMissPt[index-1]->Draw();
c_Miss_1D->cd(4);
gPad->SetLogx();
hphoton_m2[index-1]->Draw();
c_Miss_1D->cd(5);
gPad->SetLogy();
hMissMM2[index-1]->Draw();
c_Miss_1D->cd(6);
gPad->SetLogy();
hMissPt[index-1]->Draw();

TCanvas *c_flux_factor = new TCanvas("flux_factor","flux_factor",600,600);
c_flux_factor->Divide(2,2);
c_flux_factor->cd(1);
hEg->Draw();
c_flux_factor->cd(2);
hflux_factor->Draw();
c_flux_factor->cd(3);
hEgflux_factor->Draw();
c_flux_factor->cd(4);
hEg_flux_factor->Draw("colz");

// cout << "hEgflux_factor " << hEgflux_factor->Integral() << endl;

TCanvas *c_kin = new TCanvas("kin","kin",1800,600);
c_kin->Divide(9,2);
c_kin->cd(1);
gPad->SetLogy();
hs_final->Draw();
c_kin->cd(2);
gPad->SetLogy();
heta_final->Draw();
c_kin->cd(3);
gPad->SetLogy();
ht_final->Draw();
c_kin->cd(4);
gPad->SetLogy();
hQ2_final->Draw();
c_kin->cd(5);
gPad->SetLogz();
hs_t_final->Draw("colz");
c_kin->cd(6);
gPad->SetLogz();
hs_Q2_final->Draw("colz");
c_kin->cd(7);
gPad->SetLogz();
heta_t_final->Draw("colz");
c_kin->cd(8);
gPad->SetLogz();
heta_Q2_final->Draw("colz");
c_kin->cd(9);
gPad->SetLogz();
ht_Q2_final->Draw("colz");
c_kin->cd(10);
hs_final_bin->Draw();
c_kin->cd(11);
heta_final_bin->Draw();
c_kin->cd(12);
ht_final_bin->Draw();
c_kin->cd(13);
hQ2_final_bin->Draw();
c_kin->cd(14);
hs_t_final_bin->Draw("colz");
c_kin->cd(15);
hs_Q2_final_bin->Draw("colz");
c_kin->cd(16);
heta_t_final_bin->Draw("colz");
c_kin->cd(17);
heta_Q2_final_bin->Draw("colz");
c_kin->cd(18);
ht_Q2_final_bin->Draw("colz");

TCanvas *c_kin_3D = new TCanvas("kin_3D","kin_3D",1600,800);
c_kin_3D->Divide(2,2);
c_kin_3D->cd(1);
heta_t_Q2_final->Draw();
c_kin_3D->cd(3);
heta_t_Q2_final_bin->Draw("box");
c_kin_3D->cd(2);
hs_t_Q2_final->Draw();
c_kin_3D->cd(4);
hs_t_Q2_final_bin->Draw("box");

cout << "ht_Q2_etabin" << endl;
TCanvas *c_t_Q2_etabin = new TCanvas("t_Q2_etabin","t_Q2_etabin",1600,800);
c_t_Q2_etabin->Divide(Netabin/2,2);
for(int etabin=0;etabin<Netabin;etabin++){
c_t_Q2_etabin->cd(etabin+1);
gPad->SetLogz();
ht_Q2_etabin[etabin]->Draw("colz");
char title[100];
sprintf(title,"%.03f < #eta < %.03f;t (GeV^{2});Q'2 (GeV^{2})",etabin_edge[etabin],etabin_edge[etabin+1]);
ht_Q2_etabin[etabin]->SetTitle(title);
cout << ht_Q2_etabin[etabin]->Integral() << endl;
}

cout << "etabin " << endl;
for(int etabin=0;etabin<Netabin;etabin++){
for(int Q2bin=0;Q2bin<NQ2bin_etabin;Q2bin++){
for(int tbin=0;tbin<1;tbin++){  
  cout << countbin_etabin[etabin][Q2bin][tbin] <<  "\t" ;
}}
  cout << endl;
}


TCanvas *c_theta_phi_CM_etabin = new TCanvas("theta_phi_CM_etabin","theta_phi_CM_etabin",1600,800);
c_theta_phi_CM_etabin->Divide(Netabin/2,2);
for(int etabin=0;etabin<Netabin;etabin++){
c_theta_phi_CM_etabin->cd(etabin+1);
gPad->SetLogz();
htheta_phi_CM_etabin[etabin]->Draw("colz");
char title[100];
sprintf(title,"%.03f < #eta < %.03f;#phi_{CM} (rad);#theta_{CM} (rad})",etabin_edge[etabin],etabin_edge[etabin+1]);
htheta_phi_CM_etabin[etabin]->SetTitle(title);
}
	
cout << "ht_Q2_etabin_even" << endl;
TCanvas *c_t_Q2_etabin_even = new TCanvas("t_Q2_etabin_even","t_Q2_etabin_even",1600,800);
c_t_Q2_etabin_even->Divide(Netabin/2,2);
for(int etabin_even=0;etabin_even<Netabin;etabin_even++){
c_t_Q2_etabin_even->cd(etabin_even+1);
gPad->SetLogz();
ht_Q2_etabin_even[etabin_even]->Draw("colz");
// cout << ht_Q2_etabin_even[etabin_even]->Integral() << endl;
}

cout << "ht_eta_Q2bin" << endl;
TCanvas *c_t_eta_Q2bin = new TCanvas("t_eta_Q2bin","t_eta_Q2bin",1600,400);
c_t_eta_Q2bin->Divide(NQ2bin,1);
for(int Q2bin=0;Q2bin<NQ2bin;Q2bin++){
c_t_eta_Q2bin->cd(Q2bin+1);
gPad->SetLogz();
ht_eta_Q2bin[Q2bin]->Draw("colz");
char title[100];
sprintf(title,"%.01f < Q'2 < %.01f (GeV^{2});t (GeV^{2});#eta ",Q2bin_edge[Q2bin],Q2bin_edge[Q2bin+1]);
ht_eta_Q2bin[Q2bin]->SetTitle(title);
cout << ht_eta_Q2bin[Q2bin]->Integral() << endl;
}


  cout << "Q2bin " << endl;
for(int Q2bin=0;Q2bin<NQ2bin;Q2bin++){
for(int etabin=0;etabin<Netabin_Q2bin;etabin++){
for(int tbin=0;tbin<1;tbin++){  
  cout << countbin_Q2bin[Q2bin][etabin][tbin] <<  "\t" ;
}}
  cout << endl;
}

TCanvas *c_theta_phi_CM_Q2bin = new TCanvas("theta_phi_CM_Q2bin","theta_phi_CM_Q2bin",1600,400);
c_theta_phi_CM_Q2bin->Divide(NQ2bin,1);
for(int Q2bin=0;Q2bin<NQ2bin;Q2bin++){
c_theta_phi_CM_Q2bin->cd(Q2bin+1);
gPad->SetLogz();
htheta_phi_CM_Q2bin[Q2bin]->Draw("colz");
char title[100];
sprintf(title,"%.01f < Q'2 < %.01f (GeV^{2});#phi_{CM} (rad);#theta_{CM} (rad})",Q2bin_edge[Q2bin],Q2bin_edge[Q2bin+1]);
htheta_phi_CM_Q2bin[Q2bin]->SetTitle(title);
}



cout << "ht_Q2_sbin" << endl;
TCanvas *c_t_Q2_sbin = new TCanvas("t_Q2_sbin","t_Q2_sbin",1600,800);
c_t_Q2_sbin->Divide(Nsbin/2,2);
for(int sbin=0;sbin<Nsbin;sbin++){
c_t_Q2_sbin->cd(sbin+1);
gPad->SetLogz();
ht_Q2_sbin[sbin]->Draw("colz");
char title[100];
sprintf(title,"%.02f < s < %.02f;t (GeV^{2});Q'2 (GeV^{2})",sbin_edge[sbin],sbin_edge[sbin+1]);
ht_Q2_sbin[sbin]->SetTitle(title);
cout << ht_Q2_sbin[sbin]->Integral() << endl;
}

cout << "ht_Q2_sbin_even" << endl;
TCanvas *c_t_Q2_sbin_even = new TCanvas("t_Q2_sbin_even","t_Q2_sbin_even",1600,800);
c_t_Q2_sbin_even->Divide(Nsbin/2,2);
for(int sbin_even=0;sbin_even<Nsbin;sbin_even++){
c_t_Q2_sbin_even->cd(sbin_even+1);
gPad->SetLogz();
ht_Q2_sbin_even[sbin_even]->Draw("colz");
cout << ht_Q2_sbin_even[sbin_even]->Integral() << endl;
}

TCanvas *c_t_Q2_xi = new TCanvas("t_Q2_xi","t_Q2_xi",1600,800);
c_t_Q2_xi->Divide(3,1);
c_t_Q2_xi->cd(1);
gPad->SetLogz();
ht_Q2_xi[1]->SetTitle("generated");
ht_Q2_xi[1]->Draw("box");
c_t_Q2_xi->cd(2);
gPad->SetLogz();
ht_Q2_xi[3]->SetTitle("detected");
ht_Q2_xi[3]->Draw("box");
c_t_Q2_xi->cd(3);
gPad->SetLogz();
TH3F *hacc_t_Q2_xi=(TH3F*)ht_Q2_xi[3]->Clone();
hacc_t_Q2_xi->Divide(ht_Q2_xi[1]);
hacc_t_Q2_xi->SetMinimum(0);
hacc_t_Q2_xi->SetMaximum(1);
hacc_t_Q2_xi->SetNameTitle("acc_t_Q2_xi","acceptance");
hacc_t_Q2_xi->Draw("box");

TCanvas *c_Q2_xi = new TCanvas("Q2_xi","Q2_xi",1600,800);
c_Q2_xi->Divide(3,1);
c_Q2_xi->cd(1);
gPad->SetLogz();
hQ2_xi[1]->SetTitle("generated");
hQ2_xi[1]->Draw("colz");
c_Q2_xi->cd(2);
gPad->SetLogz();
hQ2_xi[3]->SetTitle("detected");
hQ2_xi[3]->Draw("colz");
c_Q2_xi->cd(3);
gPad->SetLogz();
TH2F *hacc_Q2_xi=(TH2F*)hQ2_xi[3]->Clone();
hacc_Q2_xi->Divide(hQ2_xi[1]);
hacc_Q2_xi->SetMinimum(0);
hacc_Q2_xi->SetMaximum(1);
hacc_Q2_xi->SetNameTitle("acc_Q2_xi","acceptance");
hacc_Q2_xi->Draw("colz");

TCanvas *c_t_Q2 = new TCanvas("t_Q2","t_Q2",1600,800);
c_t_Q2->Divide(3,1);
c_t_Q2->cd(1);
gPad->SetLogz();
hQ2_t[1]->SetTitle("generated");
hQ2_t[1]->Draw("colz");
c_t_Q2->cd(2);
gPad->SetLogz();
hQ2_t[3]->SetTitle("detected");
hQ2_t[3]->Draw("colz");
c_t_Q2->cd(3);
gPad->SetLogz();
TH2F *hacc_t_Q2=(TH2F*)hQ2_t[3]->Clone();
hacc_t_Q2->Divide(hQ2_t[1]);
hacc_t_Q2->SetMinimum(0);
hacc_t_Q2->SetMaximum(1);
hacc_t_Q2->SetNameTitle("acc_t_Q2","acceptance");
hacc_t_Q2->Draw("colz");

TCanvas *c_eta_xi = new TCanvas("eta_xi","eta_xi",1600,800);
heta_xi[1]->Draw("colz");


// cout << "number of entries" << endl;
// for(int k=0;k<n;k++){ cout << k << " "  << hQ2_t[k]->GetEntries() << endl;}
cout << "integral" << endl;
for(int k=0;k<n;k++){ cout << k << " "  << hQ2_t[k]->Integral() << endl;}

// cout << "total events before quasi photon cut " << hs[index-1]->Integral() << endl;
// cout << "total events after quasi photon cut " << hs[index]->Integral() << endl;

hacc_t_Q2_xi->SetDirectory(outputfile);
hacc_t_Q2->SetDirectory(outputfile);
hacc_Q2_xi->SetDirectory(outputfile);

outputfile->Write();
outputfile->Flush();

// rootapp->Run();

 return 0;


}