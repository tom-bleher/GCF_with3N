// De-excite the A-2 remnant of GCF events with Geant4's nuclear de-excitation (evaporation,
// Fermi break-up for light nuclei, photon emission), as a replacement for the evaporation
// stage of BeAGLE in PRC 105, 034001 (doi:10.1103/PhysRevC.105.034001).
//
//   deexcite IN.hepmc3 OUT.hepmc3 [--seed S] [--json OUT.json]
//
// Every final-state (status 1) nucleus with A >= 3 is taken as an excited remnant: its excitation
// energy is E* = m - M_gs(Z, A), with m its HepMC generated mass (GCF writes the A-2 mass including
// the sampled E*) and M_gs the Geant4 ground-state mass. The nucleus is broken up in its rest frame
// and the products are boosted back with its velocity, so energy and momentum are conserved
// exactly. The products are appended as status 1, after all other particles, so the order of the
// other final-state particles is unchanged. The remnant is kept as a status 3 (documentation) entry
// with the products as its daughters: npsim (DD4hep) does not pass a status 3 particle to Geant4 and
// starts its status 1 daughters at the primary vertex. As status 2 the remnant would be a Geant4
// primary with pre-assigned decay products, which Geant4 ignores for a nuclide with a lifetime
// (10Be, 10C): it transports the ion instead and the products are never simulated.
// Products are written as ground-state species (Geant4 ends the chain in the ground state or a
// long-lived isomer; isomer levels are dropped from the PDG code).
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"

#include "G4Alpha.hh"
#include "G4Deuteron.hh"
#include "G4Electron.hh"
#include "G4ExcitationHandler.hh"
#include "G4Fragment.hh"
#include "G4Gamma.hh"
#include "G4GenericIon.hh"
#include "G4He3.hh"
#include "G4IonTable.hh"
#include "G4Neutron.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4ReactionProduct.hh"
#include "G4ReactionProductVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4Triton.hh"
#include "G4Version.hh"
#include "Randomize.hh"

using namespace HepMC3;

namespace {

int ion_pdg(int z, int a) { return 1000000000 + 10000 * z + 10 * a; }

void usage() {
  std::cerr << "usage: deexcite IN.hepmc3 OUT.hepmc3 [--seed S] [--json OUT.json]\n";
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) { usage(); return 2; }
  std::string in = argv[1], out = argv[2], json = out + ".json";
  long seed = 1;
  for (int i = 3; i < argc; i++) {
    std::string a = argv[i];
    if (a == "--seed" && i + 1 < argc) seed = std::atol(argv[++i]);
    else if (a == "--json" && i + 1 < argc) json = argv[++i];
    else { usage(); return 2; }
  }

  // particle definitions the de-excitation can produce, then the handler
  G4Gamma::Definition(); G4Electron::Definition(); G4Positron::Definition();
  G4Proton::Definition(); G4Neutron::Definition(); G4Deuteron::Definition(); G4Triton::Definition();
  G4He3::Definition(); G4Alpha::Definition();
  // outside a run manager GenericIon has neither a definition ID nor a process manager, and
  // G4IonTable then refuses to create the ions the break-up produces
  G4ParticleDefinition* gion = G4GenericIon::Definition();
  if (gion->GetParticleDefinitionID() < 0) gion->SetParticleDefinitionID();
  if (!gion->GetProcessManager()) gion->SetProcessManager(new G4ProcessManager(gion));
  G4ParticleTable::GetParticleTable()->SetReadiness();
  G4Random::setTheSeed(seed);
  auto* handler = new G4ExcitationHandler();
  handler->Initialise();

  ReaderAscii reader(in);
  if (reader.failed()) { std::cerr << "deexcite: cannot read " << in << "\n"; return 3; }
  WriterAscii* writer = nullptr;

  long n_events = 0, n_remnants = 0, n_negative = 0, n_violations = 0;
  double sum_estar = 0, max_estar = 0, max_dp = 0, max_dp_rest = 0, max_negative = 0;
  std::map<int, long> products;                 // PDG -> count
  std::map<std::string, std::map<int, long>> mult;  // "n"/"p" -> multiplicity -> events

  GenEvent evt(Units::GEV, Units::MM);
  while (true) {
    reader.read_event(evt);
    if (reader.failed()) break;
    if (!writer) writer = new WriterAscii(out, reader.run_info());
    evt.set_units(Units::GEV, Units::MM);
    n_events++;
    std::vector<GenParticlePtr> remnants;
    for (auto& p : evt.particles())
      if (p->status() == 1 && p->pid() > 1000000000 && (p->pid() / 10) % 1000 >= 3) remnants.push_back(p);
    long nn = 0, np = 0;
    for (auto& r : remnants) {
      const int z = (r->pid() / 10000) % 1000, a = (r->pid() / 10) % 1000;
      const FourVector P = r->momentum();
      const double m_in = r->is_generated_mass_set() ? r->generated_mass() : P.m();
      const double m_gs = G4NucleiProperties::GetNuclearMass(a, z) / GeV;
      double estar = m_in - m_gs;
      if (estar < 0) { n_negative++; max_negative = std::max(max_negative, -estar); estar = 0; }
      const double m_star = m_gs + estar;
      n_remnants++; sum_estar += estar; max_estar = std::max(max_estar, estar);

      // break up at rest (MeV), then boost back with the remnant's velocity
      G4Fragment frag(a, z, G4LorentzVector(0, 0, 0, m_star * GeV));
      G4ReactionProductVector* prods = handler->BreakItUp(frag);
      const double ex = P.px() / P.e(), ey = P.py() / P.e(), ez = P.pz() / P.e();
      G4LorentzVector sum, rest;
      int sum_a = 0, sum_z = 0;
      auto v = std::make_shared<GenVertex>(r->production_vertex() ? r->production_vertex()->position() : FourVector());
      v->add_particle_in(r);
      std::vector<GenParticlePtr> outs;
      for (auto* pr : *prods) {
        G4LorentzVector lv(pr->GetMomentum(), pr->GetTotalEnergy());
        rest += lv;
        lv.boost(ex, ey, ez);
        sum += lv;
        const auto* def = pr->GetDefinition();
        sum_a += def->GetBaryonNumber();
        sum_z += int(std::lround(def->GetPDGCharge() / eplus));
        int pdg = def->GetPDGEncoding();
        if (def->GetParticleType() == "nucleus" && def->GetBaryonNumber() > 1)  // drop isomer level
          pdg = ion_pdg(int(std::lround(def->GetPDGCharge() / eplus)), def->GetBaryonNumber());
        auto gp = std::make_shared<GenParticle>(FourVector(lv.px() / GeV, lv.py() / GeV, lv.pz() / GeV, lv.e() / GeV), pdg, 1);
        gp->set_generated_mass(pr->GetMass() / GeV);
        outs.push_back(gp);
        products[pdg]++;
        nn += pdg == 2112; np += pdg == 2212;
        delete pr;
      }
      delete prods;
      // electrons from internal conversion carry charge -1 and no baryon number
      long n_e = 0;
      for (auto& gp : outs) n_e += gp->pid() == 11;
      if (sum_a != a || sum_z + n_e != z) n_violations++;
      // Geant4 conserves energy at the sub-keV level in the remnant rest frame; the boost (gamma ~ 44
      // at 41 GeV/u, 117 at 110 GeV/u) scales that up in the lab
      max_dp_rest = std::max({max_dp_rest, std::fabs(rest.px()), std::fabs(rest.py()), std::fabs(rest.pz()),
                              std::fabs(rest.e() - m_star * GeV)});  // MeV
      const double scale = m_star / m_in;  // only differs from 1 for clamped E* < 0
      max_dp = std::max({max_dp, std::fabs(sum.px() / GeV - P.px() * scale), std::fabs(sum.py() / GeV - P.py() * scale),
                         std::fabs(sum.pz() / GeV - P.pz() * scale), std::fabs(sum.e() / GeV - P.e() * scale)});
      r->set_status(3);  // history only; the products are the final state (see the header)
      evt.add_vertex(v);
      for (auto& gp : outs) v->add_particle_out(gp);
    }
    mult["n"][nn]++; mult["p"][np]++;
    writer->write_event(evt);
  }
  if (writer) writer->close();
  if (!n_events) { std::cerr << "deexcite: no events in " << in << "\n"; return 3; }

  std::ofstream js(json);
  js << "{\n  \"input\": \"" << in << "\",\n  \"output\": \"" << out << "\",\n  \"seed\": " << seed
     << ",\n  \"geant4\": \"" << G4Version << "\",\n  \"events\": " << n_events << ",\n  \"remnants\": " << n_remnants
     << ",\n  \"mean_Estar_GeV\": " << (n_remnants ? sum_estar / n_remnants : 0) << ",\n  \"max_Estar_GeV\": " << max_estar
     << ",\n  \"remnants_below_ground_state\": " << n_negative << ",\n  \"max_below_ground_state_GeV\": " << max_negative
     << ",\n  \"max_4momentum_imbalance_rest_frame_MeV\": " << max_dp_rest << ",\n  \"max_4momentum_imbalance_lab_GeV\": " << max_dp << ",\n  \"A_or_Z_not_conserved\": " << n_violations << ",\n  \"products_by_pdg\": {";
  bool first = true;
  for (auto& [pdg, c] : products) { js << (first ? "" : ", ") << "\"" << pdg << "\": " << c; first = false; }
  js << "},\n  \"events_by_multiplicity\": {";
  first = true;
  for (auto& [k, h] : mult) {
    js << (first ? "" : ", ") << "\"" << k << "\": {";
    bool f2 = true;
    for (auto& [m, c] : h) { js << (f2 ? "" : ", ") << "\"" << m << "\": " << c; f2 = false; }
    js << "}"; first = false;
  }
  js << "}\n}\n";
  std::cout << "deexcite: " << n_events << " events, " << n_remnants << " remnants, mean E* "
            << 1e3 * (n_remnants ? sum_estar / n_remnants : 0) << " MeV, max 4-momentum imbalance " << max_dp_rest
            << " MeV (rest frame), A/Z not conserved in " << n_violations << "; " << n_negative << " below the Geant4 ground state (max " << 1e3 * max_negative << " MeV)\n";
  if (n_violations || max_dp_rest > 0.01) {
    std::cerr << "deexcite: conservation check failed\n";
    return 4;
  }
  return 0;
}
