#!/usr/bin/env python3
"""Add the per-nucleon beam particle EICrecon needs to afterburned GCF events (#158).

    add_beam_nucleon.py IN.hepmc3 OUT.hepmc3

EICrecon takes the hadron beam from the status-4 protons and neutrons only (Beam.h), so with the nucleus as the
only hadron beam its beam-dependent outputs (InclusiveKinematics*, HadronicFinalState, the Breit frame) are empty.
As in ePIC's nuclear samples, each event gets an 'active nucleon': a status-4 proton along the (crossed, smeared)
nucleus beam with the nucleus momentum / A, as a second incoming particle of the beam vertex. It runs after abconv,
which accepts exactly two beam particles. Run inside eic-shell (HepMC3 through ROOT).

Caveat: the Roman Pot / OMD matrix reconstruction then finds a beam momentum (41 GeV/c at 5x41) and applies the
proton-beam matrix for it, which does not describe the Z/A = 1/2 beamline: ForwardRomanPot*RecParticles and
ForwardOffMRecParticles are not valid for these runs. The GCF analysis uses their hits only.
"""
import sys

import ROOT

M_P = 0.93827208816
ROOT.gSystem.Load("libHepMC3")
ROOT.gInterpreter.Declare(r"""
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include <cmath>
#include <stdexcept>
long add_beam_nucleon(const char* in, const char* out, double m_p) {
  using namespace HepMC3;
  ReaderAscii rd(in);
  WriterAscii* wr = nullptr;
  GenEvent ev;
  long n = 0;
  while (rd.read_event(ev) && !rd.failed()) {
    if (!wr) wr = new WriterAscii(out, rd.run_info());
    GenParticlePtr ion;
    int n_beam = 0;
    for (auto& p : ev.particles())
      if (p->status() == 4) { n_beam++; if (p->pid() > 1000000000) ion = p; }
    if (!ion || n_beam != 2) throw std::runtime_error("event without exactly one nucleus and one electron beam");
    const int a = (ion->pid() / 10) % 1000;
    const FourVector P = ion->momentum();
    const double px = P.px() / a, py = P.py() / a, pz = P.pz() / a;
    auto nuc = std::make_shared<GenParticle>(FourVector(px, py, pz, std::sqrt(px * px + py * py + pz * pz + m_p * m_p)), 2212, 4);
    nuc->set_generated_mass(m_p);
    if (!ion->end_vertex()) throw std::runtime_error("nucleus beam has no vertex");
    ion->end_vertex()->add_particle_in(nuc);
    wr->write_event(ev);
    n++;
  }
  if (wr) wr->close();
  return n;
}
""")
if len(sys.argv) != 3:
    sys.exit(__doc__)
n = ROOT.add_beam_nucleon(sys.argv[1], sys.argv[2], M_P)
if n <= 0:
    sys.exit(f"add_beam_nucleon: no events in {sys.argv[1]}")
print(f"add_beam_nucleon: {n} events, active-nucleon beam proton added")
