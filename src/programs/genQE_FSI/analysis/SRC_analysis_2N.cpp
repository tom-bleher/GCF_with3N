#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "QEGeneratorFSI.hh"
#include "fsi/GenieFSIHelpers.hh"
#include "nucleus/gcfNucleus.hh"
#include "constants.hh"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

eNCrossSection * myCS;
QEGeneratorFSI * myGen;
TRandom3 * myRand;
gcfNucleus * myNucleus;

/*
TTree-based 2N SRC analysis with FSI.
Generates events and saves ALL per-event kinematics to a ROOT TTree.
All cuts, histograms, and plotting are done in Python (analyze_2N.py).

Two physically distinct samples are produced, selected by wf_mode:

  wf_mode 0 (default):  2N SRC pair, p_rel ∈ [ generator default ]  — AV18 S(p)
  wf_mode 2 (SRC):      2N SRC pair, p_rel ∈ [ kF, 1.05 ] GeV/c     — AV18 S(p)
  wf_mode 4 (MF_1N):    SINGLE struck nucleon drawn from a Fermi gas (no
                        correlated pair), transported through GENIE FSI in the
                        A-1 residual. The nucleon momentum |p|<kF follows a
                        global flat n(p)=3/(4π kF^3) or a local Fermi gas (see
                        argv[8]). No gen-level recoil is emitted (rec_type=0,
                        recoil 4-vector zero): any (e,e'pp) second proton must
                        come from FSI secondaries — the correct mean-field
                        picture. Generated via QEGeneratorFSI::generate_event_MF.

The non-MF (AV18, k>kF) and MF samples are NOT mixed here; the paper-figure
analysis (natalie_paper_plots.py) uses them in separate figures.

CLI (positional; later args optional):

  argv[1] = n_events      (default 1000000)
  argv[2] = doFSI         (0/1; default 1)
  argv[3] = fsiModel      ("hN"/"hA"; default hN)
  argv[4] = output_path   (default events_2N.root)
  argv[5] = sigmaCM       (GeV/c; default 0.150; ignored for wf_mode 4)
  argv[6] = wf_mode       (0=AV18 full, 2=SRC-only, 4=MF single-nucleon)
  argv[7] = Ebeam         (GeV; default 5.01 = CLAS/Hall-B Ref [28])
  argv[8] = fg_mode       (wf_mode 4 only: "global" [default] or "local")
  argv[9] = kF            (GeV/c; default 0.25; SRC/MF p_rel split, MF FG
                           ceiling, and the nAboveKF counting threshold)
  argv[11] = vertex_rho_power (1-3; default 2. Radial weight rho^n r^2 of the
                           SRC-pair FSI vertex. Use 1 for the Wright et al.
                           paper convention: pair position follows the
                           one-body density)
  argv[10] = fsi_indep    (0/1; default 0. 1 = independent FSI: each nucleon
                           transported through its own GENIE record with an
                           undepleted A-2 remnant, instead of the default
                           shared depleting remnant in a single record)

A companion <output>.meta.txt file is written so the Python analysis can
recover wf_mode, kF, total weight, etc.
*/

int main(int argc, char **argv) {
    long long n_target_events = 1000000;
    if (argc > 1) n_target_events = std::stoll(argv[1]);

    bool doFSI = true;
    if (argc > 2) doFSI = (std::atoi(argv[2]) != 0);

    FSIModel fsiModel = kHN2018;
    if (argc > 3) {
        std::string modelStr = argv[3];
        if (modelStr == "hA" || modelStr == "HA") fsiModel = kHA2018;
    }

    std::string output_filename = "events_2N.root";
    if (argc > 4) output_filename = argv[4];

    double sigCM_val = 0.150;  // matches Wright et al. paper page 3
    if (argc > 5) sigCM_val = std::atof(argv[5]);

    int wf_mode = 0;
    if (argc > 6) wf_mode = std::atoi(argv[6]);

    double Ebeam = 5.01;  // GeV — CLAS/Hall-B (Ref [28] in Wright et al.)
    if (argc > 7) Ebeam = std::atof(argv[7]);

    // fg_mode (wf_mode 4 only): global flat Fermi gas vs local Fermi gas.
    std::string fg_mode_str = "global";
    if (argc > 8) fg_mode_str = argv[8];
    const bool useLocalFG = (fg_mode_str == "local" || fg_mode_str == "Local" || fg_mode_str == "LOCAL");

    if (wf_mode != 0 && wf_mode != 2 && wf_mode != 4) {
        std::cerr << "ERROR: wf_mode must be 0 (AV18 full), 2 (SRC-only), or 4 (MF single-nucleon); got "
                  << wf_mode << std::endl;
        return 1;
    }
    const bool isMF = (wf_mode == 4);
    const char* wf_mode_name = (wf_mode == 0) ? "AV18_full"
                             : (wf_mode == 2) ? "SRC_only"
                                              : "MF_1N";

    double kF = 0.25;   // GeV/c — Fermi momentum (SRC/MF p_rel split; MF FG ceiling)
    if (argc > 9) kF = std::atof(argv[9]);
    if (kF <= 0. || kF > 0.5) {
        std::cerr << "ERROR: kF (argv[9]) must be in (0, 0.5] GeV/c; got " << kF << std::endl;
        return 1;
    }

    bool fsi_indep = false;
    if (argc > 10) fsi_indep = (std::atoi(argv[10]) != 0);

    int vertex_rho_power = 2;
    if (argc > 11) vertex_rho_power = std::atoi(argv[11]);
    fsi::SetSRCVertexRhoPower(vertex_rho_power);

    const char* fsi_model_name = (fsiModel == kHN2018) ? "hN" : "hA";
    std::string fsi_backend_str = std::string("ENABLED (GENIE ") + fsi_model_name + " intranuke cascade)";

    myRand = new TRandom3(0);
    ffModel ffMod = kelly;
    csMethod csMeth = onshell;  // matches Wright et al. paper (Rosenbluth prescription)
    myCS = new eNCrossSection(csMeth, ffMod);

    // Carbon-12 (GCF on A=12, FSI on A-2=10 residual)
    myNucleus = new gcfNucleus(6, 6, (char*)"AV18");

    myNucleus->set_sigmaCM(sigCM_val);

    // Optional contact override for model-replication studies:
    //   GCF_CUSTOM_CONTACTS="Cpp0,Cnn0,Cpn0,Cpn1"
    // e.g. N. Wright's corrupted LocalFGM "Default" set loads
    // Cpn1=46.367, Cpn0=3.633, Cpp0=50 (=> ~67% pp/nn pairs).
    if (const char* cc = std::getenv("GCF_CUSTOM_CONTACTS")) {
        double cpp0, cnn0, cpn0, cpn1;
        if (sscanf(cc, "%lf,%lf,%lf,%lf", &cpp0, &cnn0, &cpn0, &cpn1) == 4) {
            myNucleus->setCustomValues(sigCM_val, myNucleus->get_Estar(),
                                       cpp0, cnn0, cpn0, cpn1);
            std::cout << "GCF_CUSTOM_CONTACTS active: Cpp0=" << cpp0
                      << " Cnn0=" << cnn0 << " Cpn0=" << cpn0
                      << " Cpn1=" << cpn1 << std::endl;
        } else {
            std::cerr << "ERROR: bad GCF_CUSTOM_CONTACTS format" << std::endl;
            return 1;
        }
    }

    myGen = new QEGeneratorFSI(Ebeam, myNucleus, myCS, myRand);
    myGen->EnableFSI(doFSI);
    if (doFSI) {
        myGen->SetFSIModel(fsiModel);
        myGen->SetFSITuning(250.);
        myGen->SetIndependentFSI(fsi_indep);
    }

    // ---- generator configuration per wf_mode ----
    //   mode 0 (AV18_full): default p_rel range (2N pair)
    //   mode 2 (SRC_only) : p_rel ∈ [kF, 1.05], cut at kF (2N pair)
    //   mode 4 (MF_1N)    : single struck nucleon, Fermi-gas momentum < kF,
    //                       FSI in the A-1 residual (no pair / no p_rel)
    const double p_max_sampling = 1.05;
    if (wf_mode == 2) {
        myGen->set_pRel_range(kF, p_max_sampling);
        myGen->set_pRel_cut(kF);
    } else if (isMF) {
        myGen->SetFermiMomentumFG(kF);
        myGen->SetFermiGasMode(useLocalFG ? QEGeneratorFSI::kLocalFG
                                          : QEGeneratorFSI::kGlobalFG);
    }

    std::cout << "\n========================================\n"
              << "  2N SRC Analysis — TTree output\n"
              << "========================================\n"
              << "  Nucleus:      C-12\n"
              << "  Beam energy:  " << Ebeam << " GeV\n"
              << "  FSI:          " << (doFSI ? fsi_backend_str.c_str() : "DISABLED") << "\n"
              << "  sigma_CM:     " << sigCM_val << " GeV/c\n"
              << "  kF:           " << kF << " GeV/c\n"
              << "  WF mode:      " << wf_mode << " (" << wf_mode_name << ")";
    if (isMF) std::cout << "   [single-nucleon MF, "
                        << (useLocalFG ? "local" : "global") << " Fermi gas, kF=" << kF << " GeV/c]";
    std::cout << "\n"
              << "  Events:       " << n_target_events << "\n"
              << "  Output:       " << output_filename << "\n"
              << "========================================\n\n";

    // ---- Create ROOT output file and TTree ----
    TFile *outfile = new TFile(output_filename.c_str(), "RECREATE");
    TTree *tree = new TTree("events", "2N SRC events with FSI");

    // Branch variables — scalars
    Double_t br_weight;
    Int_t br_lead_type, br_rec_type;
    Int_t br_doFSI;
    Double_t br_Q2, br_xB, br_nu, br_pmiss, br_scattering_angle;
    Int_t br_nAboveKF;
    Double_t br_rho;  // wavefunction (spectral function S(pRel))

    // 4-vectors as {px, py, pz, E}
    Double_t br_electron[4], br_lead_post[4], br_recoil_post[4], br_q[4];
    Double_t br_lead_pre[4], br_recoil_pre[4];

    // FSI event stats
    Int_t br_nSecondaries, br_nPions, br_nPiPlus, br_nPiMinus, br_nPiZero;
    Bool_t br_leadCX, br_recoilCX, br_leadElastic, br_recoilElastic;
    Bool_t br_leadAbsorbed = false, br_recAbsorbed = false;

    // Variable-length FSI secondaries
    static const Int_t kMaxSec = 50;
    Int_t br_nSec;
    Int_t br_sec_pdg[kMaxSec], br_sec_parentRole[kMaxSec], br_sec_rescatterCode[kMaxSec];
    Double_t br_sec_px[kMaxSec], br_sec_py[kMaxSec], br_sec_pz[kMaxSec], br_sec_E[kMaxSec];

    // Set up branches
    tree->Branch("weight", &br_weight, "weight/D");
    tree->Branch("lead_type", &br_lead_type, "lead_type/I");
    tree->Branch("rec_type", &br_rec_type, "rec_type/I");
    tree->Branch("doFSI", &br_doFSI, "doFSI/I");

    tree->Branch("electron", br_electron, "electron[4]/D");
    tree->Branch("lead_post", br_lead_post, "lead_post[4]/D");
    tree->Branch("recoil_post", br_recoil_post, "recoil_post[4]/D");
    tree->Branch("q", br_q, "q[4]/D");
    tree->Branch("lead_pre", br_lead_pre, "lead_pre[4]/D");
    tree->Branch("recoil_pre", br_recoil_pre, "recoil_pre[4]/D");

    tree->Branch("Q2", &br_Q2, "Q2/D");
    tree->Branch("xB", &br_xB, "xB/D");
    tree->Branch("nu", &br_nu, "nu/D");
    tree->Branch("pmiss", &br_pmiss, "pmiss/D");
    tree->Branch("scattering_angle", &br_scattering_angle, "scattering_angle/D");
    tree->Branch("nAboveKF", &br_nAboveKF, "nAboveKF/I");
    tree->Branch("rho", &br_rho, "rho/D");

    tree->Branch("nSecondaries", &br_nSecondaries, "nSecondaries/I");
    tree->Branch("nPions", &br_nPions, "nPions/I");
    tree->Branch("nPiPlus", &br_nPiPlus, "nPiPlus/I");
    tree->Branch("nPiMinus", &br_nPiMinus, "nPiMinus/I");
    tree->Branch("nPiZero", &br_nPiZero, "nPiZero/I");
    tree->Branch("leadCX", &br_leadCX, "leadCX/O");
    tree->Branch("recoilCX", &br_recoilCX, "recoilCX/O");
    tree->Branch("leadElastic", &br_leadElastic, "leadElastic/O");
    tree->Branch("recoilElastic", &br_recoilElastic, "recoilElastic/O");
    // Per-leg absorption flags (diagnostics only). Set when the ORIGINAL
    // nucleon ceased to exist in the cascade. The event is KEPT with its GCF
    // primary weight; lead_post / recoil_post then hold the leg's promoted
    // leading nucleon ejectum (zero only if no nucleon ejecta exists).
    // Analyses should select on the final state, not on these flags.
    tree->Branch("leadAbsorbed", &br_leadAbsorbed, "leadAbsorbed/O");
    tree->Branch("recAbsorbed",  &br_recAbsorbed,  "recAbsorbed/O");

    tree->Branch("nSec", &br_nSec, "nSec/I");
    tree->Branch("sec_pdg", br_sec_pdg, "sec_pdg[nSec]/I");
    tree->Branch("sec_parentRole", br_sec_parentRole, "sec_parentRole[nSec]/I");
    tree->Branch("sec_rescatterCode", br_sec_rescatterCode, "sec_rescatterCode[nSec]/I");
    tree->Branch("sec_px", br_sec_px, "sec_px[nSec]/D");
    tree->Branch("sec_py", br_sec_py, "sec_py[nSec]/D");
    tree->Branch("sec_pz", br_sec_pz, "sec_pz[nSec]/D");
    tree->Branch("sec_E", br_sec_E, "sec_E[nSec]/D");

    // Set the doFSI flag (constant for all events in this run)
    br_doFSI = doFSI ? 1 : 0;

    // Counters for summary (kF already declared at top of main for the wf_mode split)
    int event_count = 0;
    int success_count = 0;
    int events_zero_weight = 0;
    double total_weight = 0.0;

    // FSI summary accumulators
    int events_with_any_cx = 0, events_with_any_pion = 0, events_with_extra_nucleons = 0;
    double weight_with_any_cx = 0.0, weight_with_any_pion = 0.0, weight_with_extra_nucleons = 0.0;
    double weight_nucleons_above_kF_lt2 = 0.0, weight_nucleons_above_kF_eq2 = 0.0, weight_nucleons_above_kF_ge3 = 0.0;
    int events_nucleons_above_kF_lt2 = 0, events_nucleons_above_kF_eq2 = 0, events_nucleons_above_kF_ge3 = 0;

    const int progress_interval = 10000;
    std::cout << "Starting event loop: " << n_target_events << " events requested..." << std::endl;

    while (success_count < n_target_events) {
        event_count++;
        double weight;
        int lead_type, rec_type;
        TLorentzVector v_k_target, v_Lead_target, v_Rec_target, v_Am2_target;
        if (isMF)
            myGen->generate_event_MF(weight, lead_type, rec_type, v_k_target, v_Lead_target, v_Rec_target, v_Am2_target);
        else
            myGen->generate_event(weight, lead_type, rec_type, v_k_target, v_Lead_target, v_Rec_target, v_Am2_target);
        // Skip rejected events. FSI never zeroes the weight (the cascade
        // Pauli-blocks internally, and absorption keeps the primary GCF
        // weight, flagged via leadAbsorbed/recAbsorbed), so weight<=0 always
        // means a primary kinematic rejection.
        if (weight <= 0.0) { events_zero_weight++; continue; }

        // ---- Density value stored in the ROOT `rho` branch ----
        // Pair modes (0,2): AV18 spectral function S(p_rel) at the pre-FSI
        // relative momentum. MF mode (4): the Fermi-gas density n_FG(p1) at
        // the TRUE sampled struck-nucleon momentum (taken straight from the
        // generator to avoid the radiative/Coulomb reconstruction ambiguity).
        // In both cases the generator already folded the density into `weight`.
        double br_rho_final;   // density value to store in the ROOT branch
        if (isMF) {
            br_rho_final = myGen->GetLastMFRho();
        } else {
            TLorentzVector q_for_pRel(0.0, 0.0, Ebeam, Ebeam);
            q_for_pRel -= v_k_target;
            TLorentzVector pre_lead   = doFSI ? myGen->GetPreFSILead()   : v_Lead_target;
            TLorentzVector pre_recoil = doFSI ? myGen->GetPreFSIRecoil() : v_Rec_target;
            TVector3 p_lead_init(pre_lead.X() - q_for_pRel.X(),
                                 pre_lead.Y() - q_for_pRel.Y(),
                                 pre_lead.Z() - q_for_pRel.Z());
            TVector3 p_recoil_init(pre_recoil.X(), pre_recoil.Y(), pre_recoil.Z());
            double pRel_mag = 0.5 * (p_lead_init - p_recoil_init).Mag();
            br_rho_final = myNucleus->get_S(pRel_mag, lead_type, rec_type);
        }

        total_weight += weight;

        // ---- FSI event stats ----
        const auto &fsi_stats = myGen->GetLastFSIEventStats();
        br_leadCX = fsi_stats.leadChargeExchange;
        br_recoilCX = fsi_stats.recoilChargeExchange;
        br_leadElastic = fsi_stats.leadElasticLike;
        br_recoilElastic = fsi_stats.recoilElasticLike;
        br_leadAbsorbed = doFSI && myGen->IsLeadAbsorbed();
        br_recAbsorbed  = doFSI && myGen->IsRecAbsorbed();
        br_nSecondaries = fsi_stats.nSecondaries;
        br_nPions = fsi_stats.nPionsTotal;
        br_nPiPlus = fsi_stats.nPiPlus;
        br_nPiMinus = fsi_stats.nPiMinus;
        br_nPiZero = fsi_stats.nPiZero;

        // Summary counters
        if (br_leadCX || br_recoilCX) { events_with_any_cx++; weight_with_any_cx += weight; }
        if (br_nPions > 0) { events_with_any_pion++; weight_with_any_pion += weight; }

        // ---- FSI secondaries ----
        const auto &secs = myGen->GetLastFSISecondaries();
        br_nSec = std::min(static_cast<int>(secs.size()), kMaxSec);
        int nNucleonSec = 0;
        for (int i = 0; i < br_nSec; i++) {
            br_sec_pdg[i] = secs[i].pdg;
            br_sec_parentRole[i] = secs[i].parentRole;
            br_sec_rescatterCode[i] = secs[i].rescatterCode;
            br_sec_px[i] = secs[i].p4.X();
            br_sec_py[i] = secs[i].p4.Y();
            br_sec_pz[i] = secs[i].p4.Z();
            br_sec_E[i] = secs[i].p4.T();
            if (secs[i].pdg == 2212 || secs[i].pdg == 2112) nNucleonSec++;
        }
        if (nNucleonSec > 0) { events_with_extra_nucleons++; weight_with_extra_nucleons += weight; }

        // ---- Count nucleons above kF ----
        int nAboveKF = 0;
        if (v_Lead_target.Vect().Mag() > kF) nAboveKF++;
        if (v_Rec_target.Vect().Mag() > kF) nAboveKF++;
        for (const auto &sec : secs) {
            if ((sec.pdg == 2212 || sec.pdg == 2112) && sec.p4.Vect().Mag() > kF)
                nAboveKF++;
        }
        br_nAboveKF = nAboveKF;

        if (nAboveKF < 2) { weight_nucleons_above_kF_lt2 += weight; events_nucleons_above_kF_lt2++; }
        else if (nAboveKF == 2) { weight_nucleons_above_kF_eq2 += weight; events_nucleons_above_kF_eq2++; }
        else { weight_nucleons_above_kF_ge3 += weight; events_nucleons_above_kF_ge3++; }

        // ---- Kinematics ----
        TLorentzVector vbeam_target(0.0, 0.0, Ebeam, Ebeam);
        TLorentzVector q = vbeam_target - v_k_target;
        TVector3 p1 = v_Lead_target.Vect() - q.Vect();  // pmiss (reconstructed initial lead)
        double Q2 = -q.Mag2();
        double nu = vbeam_target.T() - v_k_target.T();
        double xB = (nu > 0) ? (Q2 / (2 * mN * nu)) : -1.0;
        double pmiss = p1.Mag();
        double scattering_angle = v_k_target.Vect().Theta() * 180. / M_PI;

        // ---- Fill branch variables ----
        br_weight = weight;
        br_lead_type = lead_type;
        br_rec_type = rec_type;
        br_Q2 = Q2;
        br_xB = xB;
        br_nu = nu;
        br_pmiss = pmiss;
        br_scattering_angle = scattering_angle;

        // 4-vectors: {px, py, pz, E}
        br_electron[0] = v_k_target.X(); br_electron[1] = v_k_target.Y();
        br_electron[2] = v_k_target.Z(); br_electron[3] = v_k_target.T();

        br_lead_post[0] = v_Lead_target.X(); br_lead_post[1] = v_Lead_target.Y();
        br_lead_post[2] = v_Lead_target.Z(); br_lead_post[3] = v_Lead_target.T();

        br_recoil_post[0] = v_Rec_target.X(); br_recoil_post[1] = v_Rec_target.Y();
        br_recoil_post[2] = v_Rec_target.Z(); br_recoil_post[3] = v_Rec_target.T();

        br_q[0] = q.X(); br_q[1] = q.Y();
        br_q[2] = q.Z(); br_q[3] = q.T();

        // Pre-FSI momenta
        if (doFSI) {
            TLorentzVector pre_lead = myGen->GetPreFSILead();
            TLorentzVector pre_recoil = myGen->GetPreFSIRecoil();
            br_lead_pre[0] = pre_lead.X(); br_lead_pre[1] = pre_lead.Y();
            br_lead_pre[2] = pre_lead.Z(); br_lead_pre[3] = pre_lead.T();
            br_recoil_pre[0] = pre_recoil.X(); br_recoil_pre[1] = pre_recoil.Y();
            br_recoil_pre[2] = pre_recoil.Z(); br_recoil_pre[3] = pre_recoil.T();
        } else {
            // No FSI: pre = post
            br_lead_pre[0] = br_lead_post[0]; br_lead_pre[1] = br_lead_post[1];
            br_lead_pre[2] = br_lead_post[2]; br_lead_pre[3] = br_lead_post[3];
            br_recoil_pre[0] = br_recoil_post[0]; br_recoil_pre[1] = br_recoil_post[1];
            br_recoil_pre[2] = br_recoil_post[2]; br_recoil_pre[3] = br_recoil_post[3];
        }

        // ---- Density branch: AV18 S(pRel) for pair modes, n_FG(p1) for MF ----
        // Computed above (see br_rho_final).
        br_rho = br_rho_final;

        // ---- Fill tree ----
        tree->Fill();
        success_count++;

        if (success_count % progress_interval == 0) {
            std::cout << "\r  " << success_count << " / " << n_target_events
                      << " events (" << (100.0 * success_count / n_target_events)
                      << "%)" << std::flush;
        }
    }
    std::cout << "\r  " << success_count << " / " << n_target_events
              << " events (100%)    " << std::endl;

    // ---- Write tree and close ----
    outfile->cd();
    tree->Write();
    outfile->Close();

    // ---- Companion metadata file  (<output>.meta.txt, key=value) ----
    {
        std::string meta_path = output_filename + ".meta.txt";
        std::ofstream meta(meta_path);
        if (meta.is_open()) {
            meta << "# 2N SRC generator metadata\n";
            meta << "wf_mode=" << wf_mode << "\n";
            meta << "wf_mode_name=" << wf_mode_name << "\n";
            meta << "n_events_target=" << n_target_events << "\n";
            meta << "kF=" << kF << "\n";
            meta << "fsi_indep=" << (fsi_indep ? 1 : 0) << "\n";
            if (isMF) meta << "fg_mode=" << (useLocalFG ? "local" : "global") << "\n";
            meta << "sigma_CM=" << sigCM_val << "\n";
            meta << "Ebeam=" << Ebeam << "\n";
            meta << "doFSI=" << (doFSI ? 1 : 0) << "\n";
            meta << "fsi_model=" << fsi_model_name << "\n";
            meta << "src_vertex_rho_power=" << fsi::GetSRCVertexRhoPower() << "\n";
            meta << "total_weight=" << total_weight << "\n";
            meta << "n_events_saved=" << success_count << "\n";
            meta << "n_events_attempts=" << event_count << "\n";
            meta.close();
            std::cout << "Metadata written to: " << meta_path << "\n";
        } else {
            std::cerr << "WARNING: could not write metadata file " << meta_path << "\n";
        }
    }

    // ---- Print summary ----
    double gen_eff = 100.0 * success_count / event_count;
    std::cout << "\n========================================\n"
              << "  Summary\n"
              << "========================================\n"
              << "  Events generated (total attempts): " << event_count << "\n"
              << "  Events saved (weight > 0):         " << success_count << "\n"
              << "  Events with weight = 0:            " << events_zero_weight << "\n"
              << "  Generation efficiency:             " << gen_eff << "%\n"
              << "  Total weight:                      " << total_weight << "\n"
              << "\n--- FSI summary ---\n";
    if (doFSI) {
        std::cout << "  Charge exchange:      " << events_with_any_cx
                  << " events (" << (100.0 * weight_with_any_cx / total_weight) << "% weight)\n"
                  << "  Pion production:      " << events_with_any_pion
                  << " events (" << (100.0 * weight_with_any_pion / total_weight) << "% weight)\n"
                  << "  Extra nucleon sec.:   " << events_with_extra_nucleons
                  << " events (" << (100.0 * weight_with_extra_nucleons / total_weight) << "% weight)\n";
    } else {
        std::cout << "  (FSI disabled)\n";
    }
    std::cout << "\n--- Final nucleon multiplicity above kF=" << kF << " GeV/c ---\n"
              << "  <2 above kF:   " << events_nucleons_above_kF_lt2
              << "  (wfrac " << (100.0 * weight_nucleons_above_kF_lt2 / total_weight) << "%)\n"
              << "  =2 above kF:   " << events_nucleons_above_kF_eq2
              << "  (wfrac " << (100.0 * weight_nucleons_above_kF_eq2 / total_weight) << "%)\n"
              << "  >=3 above kF:  " << events_nucleons_above_kF_ge3
              << "  (wfrac " << (100.0 * weight_nucleons_above_kF_ge3 / total_weight) << "%)\n"
              << "========================================\n"
              << "Output written to: " << output_filename << "\n";

    delete myGen;
    delete myCS;
    delete myNucleus;
    delete myRand;

    return 0;
}
