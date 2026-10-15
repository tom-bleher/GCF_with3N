#include "fsi/FSIStage.hh"
#include "constants.hh"

#include <iostream>
#include <utility>

namespace fsi {

namespace {

TLorentzVector SampleVertex(int A, int power, int nTransported, TRandom3 *rnd)
{
  if (power == kVertexAuto)
    power = (nTransported <= 1) ? 1 : (nTransported == 2 ? 2 : 3);
  if (power == 1) return SampleMFPosition(A, rnd);
  if (power == 3) return SampleSRCPosition3N(A, rnd);
  // power 2: honours SetSRCVertexRhoPower / FSI_SRC_VERTEX_RHO_POWER, like
  // the 2N generator's own vertex.
  return SampleSRCPosition(A, rnd);
}

} // namespace

void TransportEvent(const StageInput &in, TransportMode mode, FSIModel model,
                    int vertexPower, TRandom3 *rnd, StageResult &out)
{
  out = StageResult();

  // Default final state: everything unchanged. Non-nucleons are never
  // transported, so they are carried over whether or not the cascade runs.
  std::vector<FSIInputNucleon> nucleons;
  std::vector<FSIParticle> passthrough;
  int removedP = 0, removedN = 0;
  for (const FSIInputNucleon &p : in.particles) {
    FSIParticle part;
    part.pdg = p.pdg; part.p4 = p.p4; part.rescatterCode = -1;
    part.generation = 0; part.parentRole = p.parentRole;
    out.particles.push_back(part);
    if (p.pdg == pCode || p.pdg == nCode) {
      nucleons.push_back(p);
      if (p.pdg == pCode) ++removedP; else ++removedN;
    } else {
      passthrough.push_back(part);
    }
  }

  out.residualZ = in.targetZ - removedP;
  const int residualN = (in.targetA - in.targetZ) - removedN;
  out.residualA = out.residualZ + residualN;
  if (nucleons.empty() || out.residualA < 3 || out.residualZ < 0 || residualN < 0)
    return;

  if (out.residualA < 7) {
    static bool warned = false;
    if (!warned) {
      std::cerr << "Warning: GENIE intranuke applied to light residual nucleus "
                << "(A_res=" << out.residualA << "). GENIE is validated for A>6; "
                << "results for A_res<7 should be treated with caution."
                << std::endl;
      warned = true;
    }
  }

  // All nucleons leave one point, weighted by rho^n of the FULL nucleus
  // since the correlated cluster formed before knockout.
  if (in.haveVertex) {
    out.vertex = in.vertex;
  } else {
    out.vertex = SampleVertex(in.targetA, vertexPower,
                              static_cast<int>(nucleons.size()), rnd);
    out.vertexSampled = true;
  }
  for (FSIInputNucleon &n : nucleons) n.x4 = out.vertex;

  std::vector<FSIParticle> cascade;
  bool ran = false;
  if (mode == kIndependent) {
    for (const FSIInputNucleon &n : nucleons) {
      std::vector<FSIInputNucleon> one(1, n);
      ran |= ApplyGenieFSI(out.residualA, out.residualZ, one, cascade, rnd, model);
    }
  } else {
    // GENIE transports the record in insertion order, so a fixed order would
    // systematically give the first nucleon the fullest medium; no order is
    // physical, so shuffle (Fisher-Yates).
    for (size_t i = nucleons.size(); i > 1; --i) {
      const size_t j = static_cast<size_t>(rnd->Integer(static_cast<UInt_t>(i)));
      std::swap(nucleons[i - 1], nucleons[j]);
    }
    ran = ApplyGenieFSI(out.residualA, out.residualZ, nucleons, cascade, rnd, model);
  }

  if (!ran) return;
  out.ran = true;
  out.particles = cascade;
  out.particles.insert(out.particles.end(), passthrough.begin(), passthrough.end());
}

} // namespace fsi
