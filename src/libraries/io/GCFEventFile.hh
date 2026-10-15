#ifndef __GCF_EVENT_FILE_H__
#define __GCF_EVENT_FILE_H__

// The GCF event-file format shared by the pipeline stages (gcf_generate,
// gcf_fsi, ...). Full description: docs/PIPELINE.md.
//
// One TTree, default name "events", one entry per event. Units: GeV, GeV/c,
// fm; target rest frame, beam along +z.
//
//   weight                  double   event weight
//   target_A, target_Z      int      target nucleus
//   beam_E                  double   beam energy
//   e_px e_py e_pz e_E      double   scattered electron
//   nuc_pdg                 vector<int>     particles leaving the primary
//   nuc_role                vector<int>     vertex (the struck nucleon and its
//   nuc_px nuc_py nuc_pz    vector<double>  correlated partners); role 0 is
//   nuc_E                   vector<double>  the struck nucleon
//   vtx_valid               int      1 if vtx_* hold the interaction vertex
//   vtx_x vtx_y vtx_z       double   vertex inside the nucleus (fm)
//   res_px res_py res_pz res_E  double  residual nucleus (before FSI)
//
// plus, per file: TH1D "hAttempts" (1 bin, total generation attempts
// including rejected ones; hadd sums it) and a TNamed per stage with
// provenance (command line, seed, git version).
//
// gcf_fsi only needs nuc_pdg and nuc_px/py/pz/E, and target_A/target_Z
// (or the --target option); everything else is optional and passed through.

#include "TLorentzVector.h"
#include "TTree.h"
#include <string>
#include <vector>

namespace gcfio {

const char * const kTreeName = "events";
const char * const kAttemptsName = "hAttempts";

// Branch buffers for writing the stage-1 part of the format.
struct EventBuffers {
  double weight = 0.;
  int targetA = 0, targetZ = 0;
  double beamE = 0.;
  double e_px = 0., e_py = 0., e_pz = 0., e_E = 0.;
  std::vector<int> nuc_pdg, nuc_role;
  std::vector<double> nuc_px, nuc_py, nuc_pz, nuc_E;
  int vtx_valid = 0;
  double vtx_x = 0., vtx_y = 0., vtx_z = 0.;
  double res_px = 0., res_py = 0., res_pz = 0., res_E = 0.;

  void Book(TTree *t) {
    t->Branch("weight", &weight, "weight/D");
    t->Branch("target_A", &targetA, "target_A/I");
    t->Branch("target_Z", &targetZ, "target_Z/I");
    t->Branch("beam_E", &beamE, "beam_E/D");
    t->Branch("e_px", &e_px, "e_px/D");
    t->Branch("e_py", &e_py, "e_py/D");
    t->Branch("e_pz", &e_pz, "e_pz/D");
    t->Branch("e_E", &e_E, "e_E/D");
    t->Branch("nuc_pdg", &nuc_pdg);
    t->Branch("nuc_role", &nuc_role);
    t->Branch("nuc_px", &nuc_px);
    t->Branch("nuc_py", &nuc_py);
    t->Branch("nuc_pz", &nuc_pz);
    t->Branch("nuc_E", &nuc_E);
    t->Branch("vtx_valid", &vtx_valid, "vtx_valid/I");
    t->Branch("vtx_x", &vtx_x, "vtx_x/D");
    t->Branch("vtx_y", &vtx_y, "vtx_y/D");
    t->Branch("vtx_z", &vtx_z, "vtx_z/D");
    t->Branch("res_px", &res_px, "res_px/D");
    t->Branch("res_py", &res_py, "res_py/D");
    t->Branch("res_pz", &res_pz, "res_pz/D");
    t->Branch("res_E", &res_E, "res_E/D");
  }

  void ClearParticles() {
    nuc_pdg.clear(); nuc_role.clear();
    nuc_px.clear(); nuc_py.clear(); nuc_pz.clear(); nuc_E.clear();
  }

  void AddParticle(int pdg, int role, const TLorentzVector &p4) {
    nuc_pdg.push_back(pdg); nuc_role.push_back(role);
    nuc_px.push_back(p4.X()); nuc_py.push_back(p4.Y());
    nuc_pz.push_back(p4.Z()); nuc_E.push_back(p4.T());
  }

  void SetElectron(const TLorentzVector &k) {
    e_px = k.X(); e_py = k.Y(); e_pz = k.Z(); e_E = k.T();
  }

  void SetResidual(const TLorentzVector &r) {
    res_px = r.X(); res_py = r.Y(); res_pz = r.Z(); res_E = r.T();
  }

  void SetVertex(bool valid, const TLorentzVector &x4) {
    vtx_valid = valid ? 1 : 0;
    vtx_x = valid ? x4.X() : 0.; vtx_y = valid ? x4.Y() : 0.;
    vtx_z = valid ? x4.Z() : 0.;
  }
};

} // namespace gcfio

#endif
