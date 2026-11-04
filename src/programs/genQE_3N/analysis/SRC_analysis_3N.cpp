#include <iostream>
#include <cmath>
#include "QEGenerator_3N.hh"
#include "constants.hh"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

eNCrossSection * myCS;
QEGenerator_3N * myGen;
TRandom3 * myRand;

/*
TTree-based 3N SRC analysis.
Generates events and saves ALL per-event kinematics to a ROOT TTree.
All cuts, histograms, and plotting are done in Python (analyze_3N.py).

Usage: ./SRC_analysis_3N [n_events] [theta_bins] [use_CM=0] [output=events_3N.root]
*/

int n_events = 1000000;
bool use_CM = false;

int main(int argc, char **argv) {
    if (argc > 1) {
        n_events = std::stoi(argv[1]);
    }
    if (argc > 2) {
        // argv[2] was theta_bins, now unused but kept for backward compat
    }
    if (argc > 3) {
        use_CM = (std::stoi(argv[3]) == 1);
    }

    std::string output_filename = "events_3N.root";
    if (argc > 4) output_filename = argv[4];

    myRand = new TRandom3(0);
    ffModel ffMod = kelly;
    csMethod csMeth = cc1;
    myCS = new eNCrossSection(csMeth, ffMod);
    eNCrossSection *myCS_onshell = new eNCrossSection(onshell, ffMod);
    eNCrossSection *myCS_cc1 = new eNCrossSection(cc1, ffMod);
    eNCrossSection *myCS_cc2 = new eNCrossSection(cc2, ffMod);
    const double Ebeam = 6.0; // GeV
    myGen = new QEGenerator_3N(Ebeam, myCS, 1, myRand);

    std::cout << "\n========================================\n"
              << "  3N SRC Analysis — TTree output\n"
              << "========================================\n"
              << "  Nucleus:      A=" << myGen->GetTargetA()
              << ", Z=" << myGen->GetTargetZ()
              << " (residual A=" << myGen->GetTargetA() - 3
              << ", Z=" << myGen->GetTargetZ() - 2 << ")\n"
              << "  Beam energy:  " << Ebeam << " GeV\n"
              << "  Events:       " << n_events << "\n"
              << "  CM motion:    " << (use_CM ? "ON" : "OFF") << "\n"
              << "  Output:       " << output_filename << "\n"
              << "========================================\n\n";

    // ---- Create ROOT output file and TTree ----
    TFile *outfile = new TFile(output_filename.c_str(), "RECREATE");
    TTree *tree = new TTree("events", "3N SRC events"); 

    // Branch variables — scalars
    Double_t br_weight;
    Int_t br_N1_type, br_N2_type, br_N3_type;
    Double_t br_Q2, br_xB, br_nu, br_pmiss, br_scattering_angle;
    Double_t br_sigma_onshell, br_sigma_cc1, br_sigma_cc2;
    Double_t br_rho, br_delta_jacobian;

    // 4-vectors as {px, py, pz, E}
    Double_t br_electron[4], br_lead[4], br_recoil2[4], br_recoil3[4];
    Double_t br_Am3[4], br_q[4];

    // Set up branches
    tree->Branch("weight", &br_weight, "weight/D");
    tree->Branch("N1_type", &br_N1_type, "N1_type/I");
    tree->Branch("N2_type", &br_N2_type, "N2_type/I");
    tree->Branch("N3_type", &br_N3_type, "N3_type/I");

    tree->Branch("electron", br_electron, "electron[4]/D");
    tree->Branch("lead", br_lead, "lead[4]/D");
    tree->Branch("recoil2", br_recoil2, "recoil2[4]/D");
    tree->Branch("recoil3", br_recoil3, "recoil3[4]/D");
    tree->Branch("Am3", br_Am3, "Am3[4]/D");
    tree->Branch("q", br_q, "q[4]/D");

    tree->Branch("Q2", &br_Q2, "Q2/D");
    tree->Branch("xB", &br_xB, "xB/D");
    tree->Branch("nu", &br_nu, "nu/D");
    tree->Branch("pmiss", &br_pmiss, "pmiss/D");
    tree->Branch("scattering_angle", &br_scattering_angle, "scattering_angle/D");
    tree->Branch("sigma_onshell", &br_sigma_onshell, "sigma_onshell/D");
    tree->Branch("sigma_cc1", &br_sigma_cc1, "sigma_cc1/D");
    tree->Branch("sigma_cc2", &br_sigma_cc2, "sigma_cc2/D");
    tree->Branch("rho", &br_rho, "rho/D");
    tree->Branch("delta_jacobian", &br_delta_jacobian, "delta_jacobian/D");

    // Counters
    int event_count = 0;
    int success_count = 0;
    double total_weight = 0.0;
    const int progress_interval = 10000;

    std::cout << "Starting event loop: " << n_events << " events requested..." << std::endl;

    while (success_count < n_events) {
        event_count++;
        double weight, Estar = 0, rho_val = 0, delta_jac = 0;
        int N1_type, N2_type, N3_type;
        TLorentzVector v_k_target, v_Lead_target, v_2_target, v_3_target, v_Am3_target;
        myGen->generate_event(weight, N1_type, N2_type, N3_type,
                              v_k_target, v_Lead_target, v_2_target, v_3_target, v_Am3_target, use_CM, Estar, rho_val, delta_jac);
        if (weight <= 0.0) { continue; }

        total_weight += weight;

        // ---- Compute derived kinematics ----
        TLorentzVector vbeam_target(0.0, 0.0, Ebeam, Ebeam);
        TLorentzVector q = vbeam_target - v_k_target;
        TVector3 p1 = v_Lead_target.Vect() - q.Vect(); // pmiss = initial lead
        double Q2 = -q.Mag2();
        double nu = vbeam_target.T() - v_k_target.T();
        double xB = (nu > 0) ? (Q2 / (2 * mN * nu)) : -1.0;
        double pmiss = p1.Mag();
        double scattering_angle = v_k_target.Vect().Theta() * 180. / M_PI;

        // ---- Fill branch variables ----
        br_weight = weight;
        br_N1_type = N1_type;
        br_N2_type = N2_type;
        br_N3_type = N3_type;
        br_Q2 = Q2;
        br_xB = xB;
        br_nu = nu;
        br_pmiss = pmiss;
        br_scattering_angle = scattering_angle;
        bool isProton = (N1_type == 2212);
        br_sigma_onshell = myCS_onshell->sigma_eN(Ebeam, v_k_target.Vect(), v_Lead_target.Vect(), isProton);
        br_sigma_cc1 = myCS_cc1->sigma_eN(Ebeam, v_k_target.Vect(), v_Lead_target.Vect(), isProton);
        br_sigma_cc2 = myCS_cc2->sigma_eN(Ebeam, v_k_target.Vect(), v_Lead_target.Vect(), isProton);
        br_rho = rho_val;
        br_delta_jacobian = delta_jac;

        // 4-vectors: {px, py, pz, E}
        br_electron[0] = v_k_target.X(); br_electron[1] = v_k_target.Y();
        br_electron[2] = v_k_target.Z(); br_electron[3] = v_k_target.T();

        br_lead[0] = v_Lead_target.X(); br_lead[1] = v_Lead_target.Y();
        br_lead[2] = v_Lead_target.Z(); br_lead[3] = v_Lead_target.T();

        br_recoil2[0] = v_2_target.X(); br_recoil2[1] = v_2_target.Y();
        br_recoil2[2] = v_2_target.Z(); br_recoil2[3] = v_2_target.T();

        br_recoil3[0] = v_3_target.X(); br_recoil3[1] = v_3_target.Y();
        br_recoil3[2] = v_3_target.Z(); br_recoil3[3] = v_3_target.T();

        br_Am3[0] = v_Am3_target.X(); br_Am3[1] = v_Am3_target.Y();
        br_Am3[2] = v_Am3_target.Z(); br_Am3[3] = v_Am3_target.T();

        br_q[0] = q.X(); br_q[1] = q.Y();
        br_q[2] = q.Z(); br_q[3] = q.T();

        // ---- Fill tree ----
        tree->Fill();
        success_count++;

        if (success_count % progress_interval == 0) {
            std::cout << "\r  " << success_count << " / " << n_events
                      << " events (" << (100.0 * success_count / n_events) << "%)" << std::flush;
        }
    }
    std::cout << "\r  " << success_count << " / " << n_events << " events (100%)    " << std::endl;

    // ---- Write tree and close ----
    outfile->cd();
    tree->Write();
    outfile->Close();

    // ---- Print summary ----
    double gen_eff = 100.0 * success_count / event_count;
    std::cout << "\n========================================\n"
              << "  Summary\n"
              << "========================================\n"
              << "  Events generated (total attempts): " << event_count << "\n"
              << "  Events saved (weight > 0):         " << success_count << "\n"
              << "  Generation efficiency:             " << gen_eff << "%\n"
              << "  Total weight:                      " << total_weight << "\n"
              << "========================================\n"
              << "Output written to: " << output_filename << "\n";

    delete myGen;
    delete myCS;
    delete myRand;

    return 0;
}
