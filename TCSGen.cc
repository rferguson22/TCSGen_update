/*
 * File:   TCSGen.cc
 * Author: rafopar
 *
 * Created on January 11, 2020, 2:53 PM
 */

#include <TF1.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <iomanip>
#include <fstream>
#include <TRandom2.h>
#include <TTCSCrs.h>
#include <TTCSKine.h>
#include <KinFunctions.h>
#include <TLorentzVector.h>
#include "RadiativeCorrections.h"
#include "TGenPhaseSpace.h"
#include <vector>
#include "CrossSecMaster.h"
#include "Models.h"
#include "CrossSecFormulae.h"
#include <string>

#include <cstdlib>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace KinFuncs;

/*
 *
 */
int main(int argc, char **argv)
{

    // ==================================
    // ==== Reading the input config file
    // ==================================

    ifstream inpconfig("GenOptions.dat");

    map<std::string, std::string> m_Settings;
    if (inpconfig.is_open())
    {
        while (!inpconfig.eof())
        {
            std::string Key;
            std::string Val;
            inpconfig >> Key;
            inpconfig >> Val;
            m_Settings[Key] = Val;
            //cout<<setw(10)<<Key<<setw(20)<<m_Settings[Key]<<endl;
        }
    }
    else
    {
        cout << "Can not open the file GenOptions.dat" << endl;
        cout << "So can not initialize settings " << endl;
        cout << "Exiting" << endl;
        exit(1);
    }

    int Nsim;
    int n_perfile;
    double Eb;
    double t_lim;
    double Eg_min;
    double Eg_max;
    double MinvMin;
    bool isLund;
    double q2_cut;
    int seed;
    double vz_max;
    double vz_min;
    bool Rad_corr;
    double rad_cut_off_min;
    double rad_cut_off_max;
    

    for (map<std::string, std::string>::iterator it = m_Settings.begin(); it != m_Settings.end(); it++)
    {

        std::string key = (*it).first;
        std::string val = (*it).second;

        if (key.compare("Nsim") == 0)
        {
            Nsim = atoi(val.c_str());
        }
        else if (key.compare("NPerFile") == 0)
        {
            n_perfile = atoi(val.c_str());
        }
        else if (key.compare("Eb") == 0)
        {
            Eb = atof(val.c_str());
        }
        else if (key.compare("tLim") == 0)
        {
            t_lim = atof(val.c_str());
        }
        else if (key.compare("EgMin") == 0)
        {
            Eg_min = atof(val.c_str());
        }
        else if (key.compare("EgMax") == 0)
        {
            Eg_max = atof(val.c_str());
        }
        else if (key.compare("MinvMin") == 0)
        {
            MinvMin = atof(val.c_str());
        }
        else if (key.compare("Q2Cut") == 0)
        {
            q2_cut = atof(val.c_str());
        }
        else if (key.compare("LUND") == 0)
        {
            isLund = atoi(val.c_str());
        }
        else if (key.compare("Seed") == 0)
        {
            seed = atoi(val.c_str());
        }
        else if (key.compare("vzMax") == 0)
        {
            vz_max = atof(val.c_str());
        }
        else if (key.compare("vzMin") == 0)
        {
            vz_min = atof(val.c_str());
        }
        else if (key.compare("RAD_CORR") == 0)
        {
            Rad_corr = atof(val.c_str());
        }
        else if (key.compare("rad_cut_off_min") == 0)
        {
            rad_cut_off_min = atof(val.c_str());
        }
    }

    cout << "Valid targets are: ";
    for (const string& name : valid_targets){
	cout << name << " ";
    }
    cout << endl;

    string target_input;
    cout << "Enter the target type: ";
    cin  >> target_input;

    if(find(valid_targets.begin(),valid_targets.end(),target_input)==valid_targets.end()){
	cout << "Invalid target type." << endl;
	return 0;
    }

    cout << "Valid models for "<< target_input<<" are: ";
    vector<string> valid_models = get_valid_models();
    for (const string& type : valid_models){
	cout<<type<<" ";
    } 
    cout<<endl;

    string model_input;
    cout>>"Enter the model: ";
    cin>>model_input;

    string c_sec_input;
    cout>>"Enter the cross section type: ";
    cin>>c_sec_input;

    if(find(valid_c_sec.begin(),valid_c_sec.end(),c_sec_input)==valid_c_sec.end()){
        cout << "Invalid cross section type." << endl;
        return 0;
    }
    
    shared_ptr<cross_section> target=create_target(c_sec_input,model_input);


    cout << "Nsim = " << Nsim << endl;
    cout << "Eb = " << Eb << endl;
    cout << "t_lim = " << t_lim << endl;
    cout << "Eg_min = " << Eg_min << endl;
    cout << "Eg_max = " << Eg_max << endl;
    cout << "MinvMin = " << MinvMin << endl;
    cout << "q2_cut = " << q2_cut << endl;
    cout << "vz_max = " << vz_max << endl;
    cout << "vz_min = " << vz_min << endl;
    cout << "IsLund = " << isLund << endl;
    cout << "Rad_corr = " << Rad_corr << endl;
    cout << "rad_cut_off_min = " << rad_cut_off_min << endl;


    cout << "**************************************************" << endl;
    cout << "*******"
         << " RandomSeedActuallyUsed: " << seed << " *******" << endl;
    cout << "**************************************************" << endl;

    const double PI = 3.14159265358979312;
    const double radian = 57.2957795130823229;
    const double Mp = 0.9383;
    const double Me = 0.00051;
  //  const double Mtar=Mp;

    
    double  Mtar = target->mass;

    //  const double Minv_min = sqrt(Mtar*Mtar + 2*Mtar*Eg_min ) - Mtar;

    const double Minv_Egmin = sqrt(Mtar * Mtar + 2 * Mtar * Eg_min) - Mtar; // This is the maximum mass square that is accessible with a given Egmin,
                                                                    // if the User specified MinvMin is above this value, then EgMin needs to be overwritten to a Eg, that will allow MinvMin production.

    if (Minv_Egmin < MinvMin)
    {
        double EgMinOld = Eg_min;
        Eg_min = (MinvMin * MinvMin + 2 * Mtar * MinvMin) / (2 * Mtar);
        cout << "With the given Eg_min, the mass " << MinvMin << " GeV is not reachable, so the Eg min will be adjusted from " << EgMinOld << " GeV to " << Eg_min << " GeV" << endl;
    }

    const double MinvMin2 = MinvMin * MinvMin;

    TRandom2 rand;
    rand.SetSeed(seed);

    TTCSKine tcs_kin1(Mtar, Eb);
    TTCSCrs crs_lmlp;

//    TLorentzVector target(0., 0., 0., Mtar);
    TLorentzVector Lcm;

   // TLorentzVector beam(0.,0.,Eb,Eb);
   // TLorentzVector W=target+beam;

    //std::vector<double> masses1234 {0.};
   // double masses2[3]={Me,Me,Mtar};    

    //TGenPhaseSpace event;
    //event.SetDecay(W,3,masses2);
   // event.SetDecay(W,masses.size(),&masses[0]);
  
    bool write_root = !isLund;
    TFile *file_out;
    ofstream Lund_out;
    int file_number = 0;
    if (!isLund)
    {
        file_out = new TFile("tcs_gen.root", "Recreate");
    }
    else
    {
        Lund_out.open("TCSGen.dat", ofstream::out);
    }

    TH2D *h_ph_h_ph_cm1 = new TH2D("h_ph_h_ph_cm1", "", 200, 0., 360., 200, 0., 360.);
    TH2D *h_th_g_th_cm1 = new TH2D("h_th_g_th_cm1", "", 200, 0., 180., 200, 0., 180.);

    //================= Definition of Tree Variables =================
    double Eg, Minv, t, Q2, s, eta;
    double psf, crs_BH, crs_INT, crs_int;
    double psf_flux, flux_factor;
    TLorentzVector L_em, L_ep, L_prot, L_rad_1, L_rad_2;
    TLorentzVector L_gprime;

    double px_prot, py_prot, pz_prot, E_prot;
    double px_ep, py_ep, pz_ep, E_ep;
    double px_em, py_em, pz_em, E_em;
    double px_rad_em, py_rad_em, pz_rad_em, E_rad_em;
    double px_rad_ep, py_rad_ep, pz_rad_ep, E_rad_ep;
    double px_rad, py_rad, pz_rad, E_rad;
    double Theta_rad, Phi_rad, Angle_g_lep;
    double E_rad_cm, Theta_rad_cm, Phi_rad_cm, Angle_g_lep_cm;
    double Inv_Mass;

    TTree *tr1 = new TTree("tr1", "TCS MC events");
    tr1->Branch("L_em", "TLorentzVector", &L_em, 3200, 99);
    tr1->Branch("L_ep", "TLorentzVector", &L_ep, 3200, 99);
    tr1->Branch("L_prot", "TLorentzVector", &L_prot, 3200, 99);
    tr1->Branch("Eg", &Eg, "Eg/D");
    tr1->Branch("Q2", &Q2, "Q2/D");
    tr1->Branch("t", &t, "t/D");
    tr1->Branch("s", &s, "s/D");
    tr1->Branch("eta", &eta, "eta/D");
    tr1->Branch("psf", &psf, "psf/D");
    tr1->Branch("flux_factor", &flux_factor, "flux_factor/D");
    tr1->Branch("crs_BH", &crs_BH, "crs_BH/D");
    tr1->Branch("crs_INT", &crs_INT, "crs_INT/D");

    tr1->Branch("px_prot", &px_prot, "px_prot/D");
    tr1->Branch("py_prot", &py_prot, "py_prot/D");
    tr1->Branch("pz_prot", &pz_prot, "pz_prot/D");
    tr1->Branch("E_prot", &E_prot, "E_prot/D");
    tr1->Branch("px_ep", &px_ep, "px_ep/D");
    tr1->Branch("py_ep", &py_ep, "py_ep/D");
    tr1->Branch("pz_ep", &pz_ep, "pz_ep/D");
    tr1->Branch("E_ep", &E_ep, "E_ep/D");
    tr1->Branch("px_em", &px_em, "px_em/D");
    tr1->Branch("py_em", &py_em, "py_em/D");
    tr1->Branch("pz_em", &pz_em, "pz_em/D");
    tr1->Branch("E_em", &E_em, "E_em/D");
    tr1->Branch("px_rad", &px_rad, "px_rad/D");
    tr1->Branch("py_rad", &py_rad, "py_rad/D");
    tr1->Branch("pz_rad", &pz_rad, "pz_rad/D");
    tr1->Branch("px_rad_em", &px_rad_em, "px_rad_em/D");
    tr1->Branch("py_rad_em", &py_rad_em, "py_rad_em/D");
    tr1->Branch("pz_rad_em", &pz_rad_em, "pz_rad_em/D");
    tr1->Branch("px_rad_ep", &px_rad_ep, "px_rad_ep/D");
    tr1->Branch("py_rad_ep", &py_rad_ep, "py_rad_ep/D");
    tr1->Branch("pz_rad_ep", &pz_rad_ep, "pz_rad_ep/D");
    tr1->Branch("E_rad", &E_rad, "E_rad/D");
    tr1->Branch("E_rad_cm", &E_rad_cm, "E_rad_cm/D");
    tr1->Branch("Inv_Mass", &Inv_Mass, "Inv_Mass/D");

    tr1->Branch("Theta_rad", &Theta_rad, "Theta_rad/D");
    tr1->Branch("Phi_rad", &Phi_rad, "Phi_rad/D");
    tr1->Branch("Angle_g_lep", &Angle_g_lep, "Angle_g_lep/D");

    tr1->Branch("Theta_rad_cm", &Theta_rad_cm, "Theta_rad_cm/D");
    tr1->Branch("Phi_rad_cm", &Phi_rad_cm, "Phi_rad_cm/D");
    tr1->Branch("Angle_g_lep_cm", &Angle_g_lep_cm, "Angle_g_lep_cm/D");

    /////////////Set Rad Corr Parameters//////////////
    RadiativeCorrections Rad_corr_1(rad_cut_off_min);
    /////////////////////////////////////////////////

    for (int i = 0; i < Nsim; i++)
    {
        if (i % 1000 == 0)
        {
            cout.flush() << "Processed " << i << " events, approximetely " << double(100. * i / double(Nsim)) << "%\r";
        }

        if ((i + 1) % n_perfile == 0)
        {
            if (isLund)
            {
                Lund_out.close();
                file_number++;
                // Lund_out.open(Form("JPsi_gen_%d.txt", file_number), ofstream::out);
                Lund_out.open(Form("TCS_gen_%d.txt", file_number), ofstream::out);
            }
        }
        
        double psf_Eg = Eg_max - Eg_min;
        Eg = rand.Uniform(Eg_min, Eg_min + psf_Eg);
        flux_factor = N_EPA(Eb, Eg, q2_cut) + N_Brem(Eg, Eb);
        s = Mtar * Mtar + 2 * Mtar * Eg;
        double t_min = T_min(0., Mtar * Mtar, MinvMin2, Mtar * Mtar, s);
        double t_max = T_max(0., Mtar * Mtar, MinvMin2, Mtar * Mtar, s);
        double psf_t = t_min - TMath::Max(t_max, t_lim);

        if (t_min > t_lim)
        {
            t = rand.Uniform(t_min - psf_t, t_min);
            double Q2max = 2 * Mtar * Eg + t - (Eg / Mtar) * (2 * Mtar * Mtar - t - sqrt(t * t - 4 * Mtar * Mtar * t)); // Page 182 of my notebook. Derived using "Q2max = s + t - 2Mtar**2 + u_max" relation

            double psf_Q2 = Q2max - MinvMin2;

            Q2 = rand.Uniform(MinvMin2, MinvMin2 + psf_Q2);

            double u = 2 * Mtar * Mtar + Q2 - s - t;
            double th_qprime = acos((s * (t - u) - Mtar * Mtar * (Q2 - Mtar * Mtar)) / sqrt(Lambda(s, 0, Mtar * Mtar) * Lambda(s, Q2, Mtar * Mtar))); // Byukling Kayanti (4.9)
            double th_pprime = PI + th_qprime;

            double Pprime = 0.5 * sqrt(Lambda(s, Q2, Mtar * Mtar) / s); // Momentum in c.m. it is the same for q_pr and p_pr

            Lcm.SetPxPyPzE(0., 0., Eg, Mtar + Eg);
            L_prot.SetPxPyPzE(Pprime * sin(th_pprime), 0., Pprime * cos(th_pprime), sqrt(Pprime * Pprime + Mtar * Mtar));
            L_gprime.SetPxPyPzE(Pprime * sin(th_qprime), 0., Pprime * cos(th_qprime), sqrt(Pprime * Pprime + Q2));

            double psf_cos_th = 2.; // cos(th):(-1 : 1)
            double psf_phi_cm = 2 * PI;

            double cos_th = rand.Uniform(-1., -1 + psf_cos_th);
            double sin_th = sqrt(1 - cos_th * cos_th);
            double phi_cm = rand.Uniform(0., 0. + psf_phi_cm);

            double El = sqrt(Q2) / 2.; // Energy of lepton in the rest frame of qprime
            double Pl = sqrt(El * El - Me * Me);

            L_em.SetPxPyPzE(Pl * sin_th * cos(phi_cm), Pl * sin_th * sin(phi_cm), Pl * cos_th, El);
            L_ep.SetPxPyPzE(-Pl * sin_th * cos(phi_cm), -Pl * sin_th * sin(phi_cm), -Pl * cos_th, El);

            /////////////////Radiative correction///////////////////
            //double
            Rad_corr_1.Set_Inv_Mass(sqrt(Q2));
            bool in_rad_tail = (rand.Uniform(0, 1) > Rad_corr_1.Compute_cs_correction_factor(sqrt(Q2))); // randomly choose if the photon is above cut_off_min
            // the random number runs 

            if (Rad_corr && in_rad_tail)
                Rad_corr_1.Soft_Photon_Emission(L_em, L_ep, L_rad_1, L_rad_2);
            else
            {
                L_rad_1.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
                L_rad_2.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
            }
            //////////////////////////////////////////////////////////

            E_rad_cm = (L_rad_1 + L_rad_2).E();
            Theta_rad_cm = L_rad_1.Theta();
            Phi_rad_cm = L_rad_1.Phi();
            Angle_g_lep_cm = L_rad_1.Angle(L_em.Vect());

            L_em.RotateY(th_qprime); // Rotate in order to get Z axis be antiparallel to the p_prime direction in the CM frame
            L_ep.RotateY(th_qprime); // Rotate in order to get Z axis be antiparallel to the p_prime direction in the CM frame
            L_rad_1.RotateY(th_qprime);
            L_rad_2.RotateY(th_qprime);

            L_em.Boost(L_gprime.BoostVector());    // Move to the CM Frame
            L_ep.Boost(L_gprime.BoostVector());    // Move to the CM Frame
            L_rad_1.Boost(L_gprime.BoostVector()); // Move to the CM Frame
            L_rad_2.Boost(L_gprime.BoostVector()); // Move to the CM Frame

            L_em.Boost(Lcm.BoostVector());    // Move to the Lab Frame
            L_ep.Boost(Lcm.BoostVector());    // Move to the Lab Frame
            L_rad_1.Boost(Lcm.BoostVector()); // Move to the Lab Frame
            L_rad_2.Boost(Lcm.BoostVector()); // Move to the Lab Frame

            L_gprime.Boost(Lcm.BoostVector());
            L_prot.Boost(Lcm.BoostVector());

            double psf_phi_lab = 2 * PI;
            double phi_rot = rand.Uniform(0., psf_phi_lab);

            L_prot.RotateZ(phi_rot);
            L_gprime.RotateZ(phi_rot);
            L_em.RotateZ(phi_rot);
            L_ep.RotateZ(phi_rot);
            L_rad_1.RotateZ(phi_rot);
            L_rad_2.RotateZ(phi_rot);

            Theta_rad = L_rad_1.Theta();
            Phi_rad = L_rad_1.Phi();
            Angle_g_lep = L_rad_1.Angle(L_em.Vect());

            tcs_kin1.SetLemLepLp(L_em, L_ep, L_prot);

            psf = psf_t * psf_Q2 * psf_phi_lab * psf_cos_th * psf_phi_cm * psf_Eg;

            // crs_lmlp.Set_SQ2t(s, Q2, t);
            crs_BH = target->c_sec(s, Q2, t, -1, (phi_cm * TMath::RadToDeg()), (acos(cos_th) * TMath::RadToDeg())); // -1: cros section is not weighted by L/L0

            if (isnan(crs_BH))
                continue;

            eta = Q2 / (2 * (s - Mtar * Mtar) - Q2);

            if (Q2 < 9. && -t < 0.8 && eta < 0.8)
            {
                crs_INT = crs_lmlp.Eval_INT(s, Q2, t, -1., (phi_cm * TMath::RadToDeg()), (acos(cos_th) * TMath::RadToDeg()), 2.); // the last argumen "1" is the sc_D
                // crs_INT = crs_lmlp.Eval_INT( tcs_kin1.GetPhi_cm(), tcs_kin1.GetTheta_cm(), 1.); //the last argumen "1" is the sc_D
            }
            else
            {
                crs_INT = 0;
            }

            double vz = rand.Uniform(vz_min, vz_max);

            px_prot = L_prot.Px();
            py_prot = L_prot.Py();
            pz_prot = L_prot.Pz();
            E_prot = L_prot.E();
            px_ep = L_ep.Px();
            py_ep = L_ep.Py();
            pz_ep = L_ep.Pz();
            E_ep = L_ep.E();
            px_em = L_em.Px();
            py_em = L_em.Py();
            pz_em = L_em.Pz();
            E_em = L_em.E();
            px_rad = (L_rad_1 + L_rad_2).Px();
            py_rad = (L_rad_1 + L_rad_2).Py();
            pz_rad = (L_rad_1 + L_rad_2).Pz();
            E_rad = (L_rad_1 + L_rad_2).E();
            px_rad_em = L_rad_1.Px();
            py_rad_em = L_rad_1.Py();
            pz_rad_em = L_rad_1.Pz();
            E_rad_em = L_rad_1.E();
            px_rad_ep = L_rad_2.Px();
            py_rad_ep = L_rad_2.Py();
            pz_rad_ep = L_rad_2.Pz();
            E_rad_ep = L_rad_2.E();

            Inv_Mass = (L_em + L_ep).M();

            if (write_root)
            {
                tr1->Fill();
            }
            else
            {
                // Writing Header
                Lund_out << 5 << setw(5) << 1 << setw(5) << 1 << setw(5) << 0 << " " << setw(5) << "  " << psf << " " << setw(15) << 0 << setw(15) << flux_factor << setw(15) << 0 << setw(5) << 0 << setw(5) << " " << crs_BH << "\n"; // Writing Proton
                // Writing Electron
                Lund_out << 1 << setw(5) << -1 << setw(5) << 1 << setw(7) << 11 << setw(5) << 0 << setw(5) << 0 << setw(15) << px_em << setw(15) << py_em << setw(15) << pz_em;
                Lund_out << setw(15) << E_em << setw(15) << Me << setw(15) << 0. << setw(15) << 0. << setw(15) << vz << "\n";
                // Writing Positron
                Lund_out << 2 << setw(5) << 1 << setw(5) << 1 << setw(7) << -11 << setw(5) << 0 << setw(5) << 0 << setw(15) << px_ep << setw(15) << py_ep << setw(15) << pz_ep;
                Lund_out << setw(15) << E_ep << setw(15) << Me << setw(15) << 0. << setw(15) << 0. << setw(15) << vz << "\n";
                // Writing Proton
                Lund_out << 3 << setw(5) << 1 << setw(5) << 1 << setw(7) << 2212 << setw(5) << 0 << setw(5) << 0 << setw(15) << px_prot << setw(15) << py_prot << setw(15) << pz_prot;
                Lund_out << setw(15) << L_prot.E() << setw(15) << Mtar << setw(15) << 0. << setw(15) << 0. << setw(15) << vz << "\n";
                // Writing Photons
                Lund_out << 4 << setw(5) << 0 << setw(5) << 1 << setw(7) << 22 << setw(5) << 0 << setw(5) << 0 << setw(15) << px_rad_em << setw(15) << py_rad_em << setw(15) << pz_rad_em;
                Lund_out << setw(15) << E_rad_em << setw(15) << 0.0 << setw(15) << 0. << setw(15) << 0. << setw(15) << vz << "\n";
                Lund_out << 5 << setw(5) << 0 << setw(5) << 1 << setw(7) << 22 << setw(5) << 0 << setw(5) << 0 << setw(15) << px_rad_ep << setw(15) << py_rad_ep << setw(15) << pz_rad_ep;
                Lund_out << setw(15) << E_rad_ep << setw(15) << 0.0 << setw(15) << 0. << setw(15) << 0. << setw(15) << vz << "\n";
            }
        }
        else
        {
            cout << " |t_min| > |t_lim|" << endl;
            cout << " t_min =  " << t_min << "   t_lim = " << t_lim << "  Eg = " << Eg << endl;
        }
    }

    if (write_root)
    {
        tr1->Write();
        h_ph_h_ph_cm1->Write();
        h_th_g_th_cm1->Write();

        file_out->Close();
    }

    return 0;
}
