#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
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
                           ceiling, and the detected-nucleon counting threshold)
  argv[10] = fsi_indep    (0/1; default 0. 1 = independent FSI: each nucleon
                           transported through its own GENIE record with an
                           undepleted A-2 remnant, instead of the default
                           shared depleting remnant in a single record)
  argv[11] = vertex_rho_power (1-3; default 2. Radial weight rho^n r^2 of the
                           SRC-pair FSI vertex. Use 1 for the Wright et al.
                           paper convention: pair position follows the
                           one-body density)
  argv[12] = src_onset    ("fixed" [default] or "local"; wf_mode 2 only).
                           Where the correlated-pair branch begins:
                             fixed = k_rel > kF (argv[9]), a constant;
                             local = k_rel > kF(r), the local Fermi momentum
                                     at the sampled pair vertex (N. Wright's
                                     convention). "local" makes the SRC and MF
                                     samples exactly complementary -- the MF
                                     ceiling is already kF(r) -- closing the
                                     gap that "fixed" leaves between kF(r) and
                                     kF for the 78% of vertices with kF(r)<kF.
  argv[13] = src_require_p1 (0/1; default 0; wf_mode 2 only). 1 = additionally
                           require |p1| = |k_rel + P_cm/2| > kF(r) after the
                           c.m. motion is added (Wright's second test). Off by
                           default: the GCF split belongs on the RELATIVE
                           momentum, and her double test is what leaks pairs
                           into the mean-field class. Provided so the size of
                           that leak can be measured.

Environment:

  GCF_QSQ_RANGE="MIN:MAX"  Q^2 sampling range in GeV^2 (default: the
                           QEGeneratorFSI range, 1-10 GeV^2). EIC-energy runs
                           (Ebeam = rest-frame-equivalent electron energy,
                           e.g. 440 GeV for 5x41) reach Q^2 ~ 20-25 GeV^2.

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

    std::string src_onset_str = "fixed";
    if (argc > 12) src_onset_str = argv[12];
    if (src_onset_str != "fixed" && src_onset_str != "local") {
        std::cerr << "ERROR: src_onset (argv[12]) must be \"fixed\" or \"local\"; got "
                  << src_onset_str << std::endl;
        return 1;
    }
    const bool src_onset_local = (src_onset_str == "local");

    bool src_require_p1 = false;
    if (argc > 13) src_require_p1 = (std::atoi(argv[13]) != 0);

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
        myGen->SetSrcOnset(src_onset_local ? QEGeneratorFSI::kLocalKF
                                           : QEGeneratorFSI::kFixedKF);
        myGen->SetSrcRequireP1(src_require_p1);
        // Draw |p_rel| from k^2 S(k) rather than uniformly. Same expectation
        // values, far smaller weight spread: with uniform sampling S(k) spans
        // several decades and a handful of low-k events dominate every sum
        // (N_eff/N was 2.5% for a fixed k>0.25 floor and 0.11% for the local
        // kF(r) floor, which reaches ~0.02 GeV/c at large radius).
        // GCF_PREL_IS=0 reverts to uniform sampling, so that the two can be
        // compared directly at high statistics -- they must agree on every
        // WEIGHTED distribution, differing only in their error bars.
        bool prel_is = true;
        if (const char* e = std::getenv("GCF_PREL_IS")) prel_is = (std::atoi(e) != 0);
        myGen->set_pRel_importance_sampling(prel_is);
        if (!prel_is) std::cout << "  [GCF_PREL_IS=0] uniform |p_rel| sampling\n";
    } else if (isMF) {
        myGen->SetFermiMomentumFG(kF);
        myGen->SetFermiGasMode(useLocalFG ? QEGeneratorFSI::kLocalFG
                                          : QEGeneratorFSI::kGlobalFG);
    }

    // Q^2 sampling range (all wf_modes share solve_lepton_kinematics).
    if (const char* qr = std::getenv("GCF_QSQ_RANGE")) {
        double qlo = 0., qhi = 0.;
        if (std::sscanf(qr, "%lf:%lf", &qlo, &qhi) != 2 || qlo < 0. || qhi <= qlo) {
            std::cerr << "ERROR: GCF_QSQ_RANGE must be \"MIN:MAX\" in GeV^2 with "
                         "0 <= MIN < MAX; got \"" << qr << "\"" << std::endl;
            return 1;
        }
        myGen->set_QSq_range(qlo, qhi);
    }

    std::cout << "\n========================================\n"
              << "  2N SRC Analysis — TTree output\n"
              << "========================================\n"
              << "  Nucleus:      C-12\n"
              << "  Beam energy:  " << Ebeam << " GeV\n"
              << "  FSI:          " << (doFSI ? fsi_backend_str.c_str() : "DISABLED") << "\n"
              << "  sigma_CM:     " << sigCM_val << " GeV/c\n"
              << "  kF:           " << kF << " GeV/c\n"
              << "  Q2 range:     [" << myGen->get_QSqmin() << ", " << myGen->get_QSqmax() << "] GeV^2\n"
              << "  WF mode:      " << wf_mode << " (" << wf_mode_name << ")";
    if (isMF) std::cout << "   [single-nucleon MF, "
                        << (useLocalFG ? "local" : "global") << " Fermi gas, kF=" << kF << " GeV/c]";
    std::cout << "\n";
    if (wf_mode == 2) {
        std::cout << "  SRC onset:    " << src_onset_str
                  << (src_onset_local ? "  [k_rel > kF(r) at the pair vertex]"
                                      : "  [k_rel > kF, constant]")
                  << "\n"
                  << "  |p1|>kF(r):   " << (src_require_p1 ? "REQUIRED" : "not required")
                  << "\n";
    }
    std::cout
              << "  Events:       " << n_target_events << "\n"
              << "  Output:       " << output_filename << "\n"
              << "========================================\n\n";

    // ---- Create ROOT output file and TTree ----
    TFile *outfile = new TFile(output_filename.c_str(), "RECREATE");
    TTree *tree = new TTree("events", "2N SRC events with FSI");

    // Branch variables — scalars
    Double_t br_weight;
    Int_t br_gen_lead_type, br_gen_rec_type;
    Int_t br_doFSI;
    Double_t br_Q2, br_xB, br_nu, br_pmiss, br_scattering_angle;
    Int_t br_nProtonAboveKF, br_nNucleonAboveKF;
    Double_t br_rho;  // wavefunction (spectral function S(pRel))

    // 4-vectors as {px, py, pz, E}
    Double_t br_electron[4], br_q[4];
    Double_t br_lead[4], br_recoil[4];
    Double_t br_lead_pre[4], br_recoil_pre[4];

    // Post-FSI final state: flat particle list (GENIE's cascade output
    // verbatim; the generated nucleon(s) when FSI is off).
    static const Int_t kMaxPart = 100;
    Int_t br_nPart;
    Int_t br_part_pdg[kMaxPart], br_part_parentRole[kMaxPart], br_part_rescatterCode[kMaxPart];
    Int_t br_part_generation[kMaxPart];
    Double_t br_part_px[kMaxPart], br_part_py[kMaxPart], br_part_pz[kMaxPart], br_part_E[kMaxPart];

    // Set up branches.
    //
    // DETECTION-BASED (e,e'pp) EVENT BUILDING — identical rule to the 3N
    // pipeline (make_events_from_genT.py):
    //   * detected proton = final-state particle with pdg 2212 and |p| > kF
    //   * lead   = the detected proton closest in angle to q
    //   * recoil = the highest-momentum other detected proton
    //   * pmiss  = |p_lead - q| (from the DETECTED lead)
    // The paper's lead-identification cuts (theta_pq < 8 deg, p_lead > 0.75
    // |q|), the (e,e'pp) requirement nProtonAboveKF == 2, and all other
    // kinematic cuts are applied downstream in the figure scripts.
    tree->Branch("weight", &br_weight, "weight/D");
    tree->Branch("gen_lead_type", &br_gen_lead_type, "gen_lead_type/I");
    tree->Branch("gen_rec_type", &br_gen_rec_type, "gen_rec_type/I");
    tree->Branch("doFSI", &br_doFSI, "doFSI/I");

    tree->Branch("electron", br_electron, "electron[4]/D");
    tree->Branch("q", br_q, "q[4]/D");
    tree->Branch("lead", br_lead, "lead[4]/D");
    tree->Branch("recoil", br_recoil, "recoil[4]/D");
    tree->Branch("lead_pre", br_lead_pre, "lead_pre[4]/D");
    tree->Branch("recoil_pre", br_recoil_pre, "recoil_pre[4]/D");

    tree->Branch("Q2", &br_Q2, "Q2/D");
    tree->Branch("xB", &br_xB, "xB/D");
    tree->Branch("nu", &br_nu, "nu/D");
    tree->Branch("pmiss", &br_pmiss, "pmiss/D");
    tree->Branch("scattering_angle", &br_scattering_angle, "scattering_angle/D");
    tree->Branch("nProtonAboveKF", &br_nProtonAboveKF, "nProtonAboveKF/I");
    tree->Branch("nNucleonAboveKF", &br_nNucleonAboveKF, "nNucleonAboveKF/I");
    tree->Branch("rho", &br_rho, "rho/D");

    tree->Branch("nPart", &br_nPart, "nPart/I");
    tree->Branch("part_pdg", br_part_pdg, "part_pdg[nPart]/I");
    tree->Branch("part_parentRole", br_part_parentRole, "part_parentRole[nPart]/I");
    tree->Branch("part_rescatterCode", br_part_rescatterCode, "part_rescatterCode[nPart]/I");
    // 0 = this particle IS the transported input nucleon, unchanged by the
    // cascade; >0 = an FSI product. parentRole cannot express this: a
    // knocked-out medium nucleon inherits the role of whatever hit it.
    tree->Branch("part_generation", br_part_generation, "part_generation[nPart]/I");
    tree->Branch("part_px", br_part_px, "part_px[nPart]/D");
    tree->Branch("part_py", br_part_py, "part_py[nPart]/D");
    tree->Branch("part_pz", br_part_pz, "part_pz[nPart]/D");
    tree->Branch("part_E", br_part_E, "part_E[nPart]/D");

    // Set the doFSI flag (constant for all events in this run)
    br_doFSI = doFSI ? 1 : 0;

    // Counters for summary (kF already declared at top of main for the wf_mode split)
    int event_count = 0;
    int success_count = 0;
    int events_zero_weight = 0;
    double total_weight = 0.0;

    // Summary accumulators: detected-proton multiplicity ((e,e'pp) = exactly 2)
    double weight_protons_above_kF_lt2 = 0.0, weight_protons_above_kF_eq2 = 0.0, weight_protons_above_kF_ge3 = 0.0;
    int events_protons_above_kF_lt2 = 0, events_protons_above_kF_eq2 = 0, events_protons_above_kF_ge3 = 0;

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
            // v_Lead_target / v_Rec_target are the GENERATED (pre-FSI)
            // nucleons: FSI no longer modifies the out-parameters.
            TVector3 p_lead_init(v_Lead_target.X() - q_for_pRel.X(),
                                 v_Lead_target.Y() - q_for_pRel.Y(),
                                 v_Lead_target.Z() - q_for_pRel.Z());
            TVector3 p_recoil_init(v_Rec_target.X(), v_Rec_target.Y(), v_Rec_target.Z());
            double pRel_mag = 0.5 * (p_lead_init - p_recoil_init).Mag();
            br_rho_final = myNucleus->get_S(pRel_mag, lead_type, rec_type);
        }

        total_weight += weight;

        // ---- Post-FSI final state: flat particle list ----
        const auto &parts = myGen->GetFinalStateParticles();
        br_nPart = std::min(static_cast<int>(parts.size()), kMaxPart);
        for (int i = 0; i < br_nPart; i++) {
            br_part_pdg[i] = parts[i].pdg;
            br_part_parentRole[i] = parts[i].parentRole;
            br_part_rescatterCode[i] = parts[i].rescatterCode;
            br_part_generation[i] = parts[i].generation;
            br_part_px[i] = parts[i].p4.X();
            br_part_py[i] = parts[i].p4.Y();
            br_part_pz[i] = parts[i].p4.Z();
            br_part_E[i] = parts[i].p4.T();
        }

        // ---- Electron kinematics ----
        TLorentzVector vbeam_target(0.0, 0.0, Ebeam, Ebeam);
        TLorentzVector q = vbeam_target - v_k_target;
        double Q2 = -q.Mag2();
        double nu = vbeam_target.T() - v_k_target.T();
        double xB = (nu > 0) ? (Q2 / (2 * mN * nu)) : -1.0;
        double scattering_angle = v_k_target.Vect().Theta() * 180. / M_PI;

        // ---- Detection-based (e,e'pp) event building ----
        // Detected proton = final-state proton with |p| > kF. Lead = the
        // detected proton closest in angle to q; recoil = the
        // highest-momentum other detected proton. Identical rule to the 3N
        // converter (make_events_from_genT.py).
        int nProtonAboveKF = 0, nNucleonAboveKF = 0;
        int lead_idx = -1, rec_idx = -1;
        double best_angle = 1.e9, best_rec_p = -1.;
        for (size_t i = 0; i < parts.size(); ++i) {
            const bool isP = (parts[i].pdg == 2212);
            const bool isN = (parts[i].pdg == 2112);
            if (!(isP || isN)) continue;
            if (parts[i].p4.Vect().Mag() <= kF) continue;
            nNucleonAboveKF++;
            if (!isP) continue;
            nProtonAboveKF++;
            const double ang = parts[i].p4.Vect().Angle(q.Vect());
            if (ang < best_angle) { best_angle = ang; lead_idx = (int)i; }
        }
        for (size_t i = 0; i < parts.size(); ++i) {
            if ((int)i == lead_idx) continue;
            if (parts[i].pdg != 2212) continue;
            if (parts[i].p4.Vect().Mag() <= kF) continue;
            if (parts[i].p4.Vect().Mag() > best_rec_p) {
                best_rec_p = parts[i].p4.Vect().Mag();
                rec_idx = (int)i;
            }
        }
        const TLorentzVector v_lead_det = (lead_idx >= 0) ? parts[lead_idx].p4
                                                          : TLorentzVector();
        const TLorentzVector v_rec_det  = (rec_idx >= 0)  ? parts[rec_idx].p4
                                                          : TLorentzVector();
        br_nProtonAboveKF = nProtonAboveKF;
        br_nNucleonAboveKF = nNucleonAboveKF;

        if (nProtonAboveKF < 2) { weight_protons_above_kF_lt2 += weight; events_protons_above_kF_lt2++; }
        else if (nProtonAboveKF == 2) { weight_protons_above_kF_eq2 += weight; events_protons_above_kF_eq2++; }
        else { weight_protons_above_kF_ge3 += weight; events_protons_above_kF_ge3++; }

        // pmiss from the DETECTED lead (measurement-faithful).
        TVector3 p1 = v_lead_det.Vect() - q.Vect();
        double pmiss = (lead_idx >= 0) ? p1.Mag() : -1.0;

        // ---- Fill branch variables ----
        br_weight = weight;
        br_gen_lead_type = lead_type;
        br_gen_rec_type = rec_type;
        br_Q2 = Q2;
        br_xB = xB;
        br_nu = nu;
        br_pmiss = pmiss;
        br_scattering_angle = scattering_angle;

        // 4-vectors: {px, py, pz, E}
        br_electron[0] = v_k_target.X(); br_electron[1] = v_k_target.Y();
        br_electron[2] = v_k_target.Z(); br_electron[3] = v_k_target.T();

        br_q[0] = q.X(); br_q[1] = q.Y();
        br_q[2] = q.Z(); br_q[3] = q.T();

        br_lead[0] = v_lead_det.X(); br_lead[1] = v_lead_det.Y();
        br_lead[2] = v_lead_det.Z(); br_lead[3] = v_lead_det.T();

        br_recoil[0] = v_rec_det.X(); br_recoil[1] = v_rec_det.Y();
        br_recoil[2] = v_rec_det.Z(); br_recoil[3] = v_rec_det.T();

        // Generated (pre-FSI) nucleons: the generate_event out-parameters are
        // unmodified by FSI in the flat-list scheme.
        br_lead_pre[0] = v_Lead_target.X(); br_lead_pre[1] = v_Lead_target.Y();
        br_lead_pre[2] = v_Lead_target.Z(); br_lead_pre[3] = v_Lead_target.T();
        br_recoil_pre[0] = v_Rec_target.X(); br_recoil_pre[1] = v_Rec_target.Y();
        br_recoil_pre[2] = v_Rec_target.Z(); br_recoil_pre[3] = v_Rec_target.T();

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
            if (wf_mode == 2) {
                meta << "src_onset=" << src_onset_str << "\n";
                meta << "src_require_p1=" << (src_require_p1 ? 1 : 0) << "\n";
            }
            meta << "sigma_CM=" << sigCM_val << "\n";
            meta << "Ebeam=" << Ebeam << "\n";
            meta << "QSq_min=" << myGen->get_QSqmin() << "\n";
            meta << "QSq_max=" << myGen->get_QSqmax() << "\n";
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
              << "\n--- Detected-proton multiplicity above kF=" << kF << " GeV/c ---\n"
              << "  <2 protons:   " << events_protons_above_kF_lt2
              << "  (wfrac " << (100.0 * weight_protons_above_kF_lt2 / total_weight) << "%)\n"
              << "  =2 protons:   " << events_protons_above_kF_eq2
              << "  (wfrac " << (100.0 * weight_protons_above_kF_eq2 / total_weight) << "%)  [(e,e'pp)]\n"
              << "  >=3 protons:  " << events_protons_above_kF_ge3
              << "  (wfrac " << (100.0 * weight_protons_above_kF_ge3 / total_weight) << "%)\n"
              << "========================================\n"
              << "Output written to: " << output_filename << "\n";

    delete myGen;
    delete myCS;
    delete myNucleus;
    delete myRand;

    return 0;
}
