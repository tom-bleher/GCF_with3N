// gcf_fsi -- stage 2 of the GCF pipeline: transports the particles leaving
// the primary vertex of every event in an event file through GENIE's
// intranuclear cascade (hN or hA), and writes the file back with the
// final-state particle list added. The input can come from gcf_generate or
// from any other generator that writes the GCF event-file format
// (docs/PIPELINE.md); all of its branches are kept.

#include "fsi/FSIStage.hh"
#include "fsi/GenieFSIHelpers.hh"
#include "io/GCFEventFile.hh"
#include "io/Provenance.hh"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TFile.h"
#include "TKey.h"
#include "TNamed.h"
#include "TStopwatch.h"
#include "TTree.h"

#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;
using ROOT::RVec;

namespace {

void Usage()
{
  cerr <<
    "Usage: gcf_fsi --input FILE --output FILE [options]\n"
    "\n"
    "Adds GENIE intranuclear-cascade FSI to every event of an event file and\n"
    "writes all input branches plus the final-state particle list (part_*).\n"
    "\n"
    "  -i, --input FILE         input event file (GCF event-file format)\n"
    "  -o, --output FILE        output file (must differ from the input)\n"
    "  -t, --tree NAME          tree name (default events)\n"
    "  -f, --fsi-model hN|hA    GENIE cascade model (default hN)\n"
    "      --independent        transport each nucleon in its own GENIE record\n"
    "                           with an undepleted remnant (default: shared)\n"
    "  -S, --seed S             random seed, positive (default: drawn and printed)\n"
    "      --target A:Z         target for files without target_A/target_Z\n"
    "      --vertex-power auto|1|2|3\n"
    "                           rho^n weighting of the vertex for events without\n"
    "                           one (default auto: n = 1, 2, 3 for 1, 2, >=3 nucleons)\n"
    "      --keep REGEX         input branches to copy (default: all)\n"
    "      --overwrite          replace FSI branches already in the input\n"
    "  -n, --max-events N       process only the first N events\n"
    "  -h, --help               this message\n";
}

// Final-state buffers of the current event, shared by all output columns:
// the cascade runs once per event (in the __gcf_fsi column) and the output
// columns just read these.
struct Runner {
  fsi::TransportMode mode = fsi::kSharedRemnant;
  FSIModel model = kHN2018;
  int vertexPower = fsi::kVertexAuto;
  TRandom3 rnd;
  long long total = 0, done = 0, nRan = 0, nVertexSampled = 0, progressEvery = 1;
  fsi::StageResult res;
  vector<int> pdg, role, rescatter, generation;
  vector<double> px, py, pz, E;

  bool Run(int A, int Z, const RVec<int> &ipdg, const RVec<int> &irole,
           const RVec<double> &ipx, const RVec<double> &ipy,
           const RVec<double> &ipz, const RVec<double> &iE,
           int vtxValid, double vx, double vy, double vz)
  {
    const size_t n = ipdg.size();
    if (ipx.size() != n || ipy.size() != n || ipz.size() != n || iE.size() != n)
      throw runtime_error("nuc_pdg and nuc_px/py/pz/E have different lengths in an event");
    fsi::StageInput in;
    in.targetA = A;
    in.targetZ = Z;
    for (size_t i = 0; i < n; ++i) {
      FSIInputNucleon p;
      p.pdg = ipdg[i];
      p.p4.SetPxPyPzE(ipx[i], ipy[i], ipz[i], iE[i]);
      p.parentRole = (irole.size() == n) ? irole[i] : static_cast<int>(i);
      in.particles.push_back(p);
    }
    in.haveVertex = (vtxValid != 0);
    in.vertex.SetXYZT(vx, vy, vz, 0.);
    fsi::TransportEvent(in, mode, model, vertexPower, &rnd, res);

    pdg.clear(); role.clear(); rescatter.clear(); generation.clear();
    px.clear(); py.clear(); pz.clear(); E.clear();
    for (const FSIParticle &p : res.particles) {
      pdg.push_back(p.pdg); role.push_back(p.parentRole);
      rescatter.push_back(p.rescatterCode); generation.push_back(p.generation);
      px.push_back(p.p4.X()); py.push_back(p.p4.Y());
      pz.push_back(p.p4.Z()); E.push_back(p.p4.T());
    }
    if (res.ran) ++nRan;
    if (res.vertexSampled) ++nVertexSampled;
    if (++done % progressEvery == 0 || done == total)
      cout << "Progress: " << done << "/" << total << endl;
    return res.ran;
  }
};

} // namespace

int main(int argc, char **argv)
{
  string input, output, treeName = gcfio::kTreeName, modelName = "hN",
         vertexPowerName = "auto", keep = ".*";
  unsigned long seed = 0;
  long long maxEvents = 0;
  int defaultA = -1, defaultZ = -1;
  bool independent = false, overwrite = false;

  enum { kIndependent = 1000, kTarget, kVertexPower, kKeep, kOverwrite };
  static const option longOpts[] = {
    {"input", required_argument, nullptr, 'i'},
    {"output", required_argument, nullptr, 'o'},
    {"tree", required_argument, nullptr, 't'},
    {"fsi-model", required_argument, nullptr, 'f'},
    {"seed", required_argument, nullptr, 'S'},
    {"max-events", required_argument, nullptr, 'n'},
    {"help", no_argument, nullptr, 'h'},
    {"independent", no_argument, nullptr, kIndependent},
    {"target", required_argument, nullptr, kTarget},
    {"vertex-power", required_argument, nullptr, kVertexPower},
    {"keep", required_argument, nullptr, kKeep},
    {"overwrite", no_argument, nullptr, kOverwrite},
    {nullptr, 0, nullptr, 0}};

  int c;
  while ((c = getopt_long(argc, argv, "i:o:t:f:S:n:h", longOpts, nullptr)) != -1) {
    switch (c) {
      case 'i': input = optarg; break;
      case 'o': output = optarg; break;
      case 't': treeName = optarg; break;
      case 'f': modelName = optarg; break;
      case 'S': seed = strtoul(optarg, nullptr, 10); break;
      case 'n': maxEvents = atoll(optarg); break;
      case 'h': Usage(); return 0;
      case kIndependent: independent = true; break;
      case kTarget:
        if (sscanf(optarg, "%d:%d", &defaultA, &defaultZ) != 2 ||
            defaultA < 1 || defaultZ < 0 || defaultZ > defaultA) {
          cerr << "Invalid --target (expect A:Z, e.g. 12:6).\n"; return 1;
        }
        break;
      case kVertexPower: vertexPowerName = optarg; break;
      case kKeep: keep = optarg; break;
      case kOverwrite: overwrite = true; break;
      default: Usage(); return 1;
    }
  }

  Runner runner;
  bool ok = true;
  auto fail = [&ok](const string &msg) { cerr << "Error: " << msg << "\n"; ok = false; };
  if (input.empty() || output.empty()) fail("--input and --output are required");
  if (!input.empty() && input == output) fail("--output must differ from --input");
  if (modelName == "hN") runner.model = kHN2018;
  else if (modelName == "hA") runner.model = kHA2018;
  else fail("--fsi-model must be hN or hA");
  if (vertexPowerName == "auto") runner.vertexPower = fsi::kVertexAuto;
  else if (vertexPowerName == "1" || vertexPowerName == "2" || vertexPowerName == "3")
    runner.vertexPower = atoi(vertexPowerName.c_str());
  else fail("--vertex-power must be auto, 1, 2 or 3");
  regex keepRe;
  try { keepRe = regex(keep); } catch (const regex_error &) { fail("--keep is not a valid regex"); }
  if (!getenv("GENIE")) fail("GENIE is not set; source setup_env.sh (or build/gcf_env.sh) first");
  if (!ok) { cerr << "(gcf_fsi --help for usage)\n"; return 1; }
  runner.mode = independent ? fsi::kIndependent : fsi::kSharedRemnant;
  if (runner.vertexPower == 2) fsi::SetSRCVertexRhoPower(2);
  fsi::ResolveGenieXMLPath();

  if (seed == 0) { TRandom3 seeder(0); seed = 1 + seeder.Integer(0xfffffffe); }
  cout << "Random seed: " << seed << "\n";
  runner.rnd.SetSeed(seed);

  // ---- Input -------------------------------------------------------------
  // The events may be a TTree or an RNTuple (what uproot >= 5.6 writes by
  // default); RDataFrame reads both.
  long long nEntries = -1;
  {
    TFile in(input.c_str(), "READ");
    if (in.IsZombie()) { cerr << "Error: cannot open " << input << "\n"; return 1; }
    TKey *key = in.GetKey(treeName.c_str());
    if (!key) { cerr << "Error: no '" << treeName << "' in " << input << "\n"; return 1; }
    const string cls = key->GetClassName();
    if (cls == "TTree") {
      nEntries = in.Get<TTree>(treeName.c_str())->GetEntries();
    } else if (cls.find("RNTuple") == string::npos) {
      cerr << "Error: '" << treeName << "' in " << input << " is a " << cls
           << ", not a TTree or RNTuple\n";
      return 1;
    }
  }

  ROOT::RDataFrame df(treeName, input);
  if (nEntries < 0) nEntries = *df.Count();
  runner.total = (maxEvents > 0 && maxEvents < nEntries) ? maxEvents : nEntries;
  runner.progressEvery = max(1LL, runner.total / 100);
  const vector<string> inputColumns = df.GetColumnNames();
  auto has = [&df](const string &col) { return df.HasColumn(col); };
  for (const char *col : {"nuc_pdg", "nuc_px", "nuc_py", "nuc_pz", "nuc_E"})
    if (!has(col)) { cerr << "Error: input has no '" << col << "' branch (see docs/PIPELINE.md)\n"; return 1; }
  if ((!has("target_A") || !has("target_Z")) && defaultA < 0) {
    cerr << "Error: input has no target_A/target_Z branches; give the target with --target A:Z\n";
    return 1;
  }

  ROOT::RDF::RNode node = df;
  if (maxEvents > 0) node = node.Range(maxEvents);

  // Normalize the inputs to fixed types (std::vector or C-array branches, any
  // numeric element type).
  auto vecOf = [](const string &type, const string &col) {
    return "ROOT::VecOps::RVec<" + type + ">(" + col + ".begin(), " + col + ".end())";
  };
  node = node.Define("__gcf_A", has("target_A") ? "static_cast<int>(target_A)" : to_string(defaultA));
  node = node.Define("__gcf_Z", has("target_Z") ? "static_cast<int>(target_Z)" : to_string(defaultZ));
  node = node.Define("__gcf_pdg", vecOf("int", "nuc_pdg"));
  node = node.Define("__gcf_px", vecOf("double", "nuc_px"));
  node = node.Define("__gcf_py", vecOf("double", "nuc_py"));
  node = node.Define("__gcf_pz", vecOf("double", "nuc_pz"));
  node = node.Define("__gcf_E", vecOf("double", "nuc_E"));
  if (has("nuc_role")) node = node.Define("__gcf_role", vecOf("int", "nuc_role"));
  else node = node.Define("__gcf_role", "ROOT::VecOps::RVec<int>()");  // roles = indices
  const bool haveVertex = has("vtx_x") && has("vtx_y") && has("vtx_z");
  node = node.Define("__gcf_vvalid", !haveVertex ? "0"
                     : (has("vtx_valid") ? "static_cast<int>(vtx_valid)" : "1"));
  node = node.Define("__gcf_vx", haveVertex ? "static_cast<double>(vtx_x)" : "0.");
  node = node.Define("__gcf_vy", haveVertex ? "static_cast<double>(vtx_y)" : "0.");
  node = node.Define("__gcf_vz", haveVertex ? "static_cast<double>(vtx_z)" : "0.");

  node = node.Define("__gcf_fsi",
      [&runner](int A, int Z, const RVec<int> &pdg, const RVec<int> &role,
                const RVec<double> &px, const RVec<double> &py,
                const RVec<double> &pz, const RVec<double> &E,
                int vv, double vx, double vy, double vz) {
        return runner.Run(A, Z, pdg, role, px, py, pz, E, vv, vx, vy, vz);
      },
      {"__gcf_A", "__gcf_Z", "__gcf_pdg", "__gcf_role", "__gcf_px", "__gcf_py",
       "__gcf_pz", "__gcf_E", "__gcf_vvalid", "__gcf_vx", "__gcf_vy", "__gcf_vz"});

  // ---- Output columns ----------------------------------------------------
  vector<string> outputColumns;
  bool clash = false;
  auto add = [&](const string &name, auto &&fn) {
    if (has(name)) {
      if (!overwrite) {
        cerr << "Error: input already has '" << name
             << "' (already went through FSI?); use --overwrite to replace it\n";
        clash = true;
        return;
      }
      node = node.Redefine(name, fn, {"__gcf_fsi"});
    } else {
      node = node.Define(name, fn, {"__gcf_fsi"});
    }
    outputColumns.push_back(name);
  };
  Runner &r = runner;
  add("fsi_ran", [&r](bool) { return static_cast<int>(r.res.ran); });
  add("fsi_vtx_sampled", [&r](bool) { return static_cast<int>(r.res.vertexSampled); });
  add("fsi_vtx_x", [&r](bool) { return r.res.vertex.X(); });
  add("fsi_vtx_y", [&r](bool) { return r.res.vertex.Y(); });
  add("fsi_vtx_z", [&r](bool) { return r.res.vertex.Z(); });
  add("fsi_res_A", [&r](bool) { return r.res.residualA; });
  add("fsi_res_Z", [&r](bool) { return r.res.residualZ; });
  add("part_n", [&r](bool) { return static_cast<int>(r.pdg.size()); });
  add("part_pdg", [&r](bool) { return r.pdg; });
  add("part_parentRole", [&r](bool) { return r.role; });
  add("part_rescatterCode", [&r](bool) { return r.rescatter; });
  add("part_generation", [&r](bool) { return r.generation; });
  add("part_px", [&r](bool) { return r.px; });
  add("part_py", [&r](bool) { return r.py; });
  add("part_pz", [&r](bool) { return r.pz; });
  add("part_E", [&r](bool) { return r.E; });
  if (clash) return 1;

  vector<string> columns;
  set<string> seen(outputColumns.begin(), outputColumns.end());
  for (const string &col : inputColumns) {
    if (seen.count(col) || !regex_match(col, keepRe)) continue;
    seen.insert(col);
    columns.push_back(col);
  }
  columns.insert(columns.end(), outputColumns.begin(), outputColumns.end());

  // ---- Run ---------------------------------------------------------------
  TStopwatch timer;
  ROOT::RDF::RSnapshotOptions opts;
  opts.fMode = "RECREATE";
  try {
    node.Snapshot(treeName, output, columns, opts);
  } catch (const exception &e) {
    cerr << "Error: " << e.what() << "\n";
    return 1;
  }
  timer.Stop();

  // ---- Carry over the other objects (hAttempts, provenance of earlier
  // stages, ...) and add our own provenance --------------------------------
  {
    TFile in(input.c_str(), "READ");
    TFile out(output.c_str(), "UPDATE");
    set<string> copied;
    TIter next(in.GetListOfKeys());
    while (TKey *key = static_cast<TKey *>(next())) {
      const string name = key->GetName();
      if (name == treeName || name == "gcf_fsi" || copied.count(name)) continue;
      copied.insert(name);   // keys come highest cycle first
      if (string(key->GetClassName()) == "TTree") {
        cerr << "Note: extra tree '" << name << "' in the input is not copied\n";
        continue;
      }
      TObject *obj = key->ReadObj();
      out.cd();
      obj->Write(name.c_str());
      delete obj;
    }
    ostringstream prov;
    prov << gcfio::Provenance("gcf_fsi", argc, argv, seed)
         << "input: " << input << "\n"
         << "fsi_model: " << modelName << "\n"
         << "transport: " << (independent ? "independent" : "shared") << "\n"
         << "vertex_power: " << vertexPowerName << "\n"
         << "events: " << runner.done << ", cascade ran: " << runner.nRan
         << ", vertex sampled: " << runner.nVertexSampled << "\n";
    out.cd();
    TNamed("gcf_fsi", prov.str().c_str()).Write("gcf_fsi", TObject::kOverwrite);
  }

  cout << "FSI (" << modelName << ", " << (independent ? "independent" : "shared")
       << ") applied to " << runner.done << " events in " << timer.RealTime() << " s: "
       << "cascade ran for " << runner.nRan << ", skipped for "
       << runner.done - runner.nRan << " (residual A < 3 or nothing to transport); "
       << runner.nVertexSampled << " vertices sampled.\nWrote " << output << "\n";
  return 0;
}
