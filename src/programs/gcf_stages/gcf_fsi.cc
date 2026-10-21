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

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <map>
#include <sys/stat.h>
#include <unistd.h>
#include <libxml/parser.h>
#include <iostream>
#include <regex>
#include <set>
#include <cerrno>
#include <cstring>
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
    "      --channel NAME=F     scale the rate of an interaction channel of the cascade\n"
    "                           by F (0 = switch it off); repeatable. --list-channels\n"
    "                           shows the names and which model each applies to\n"
    "      --list-channels      print the channels and exit\n"
    "      --keep REGEX         input branches to copy (default: all)\n"
    "      --overwrite          replace FSI branches already in the input\n"
    "  -n, --max-events N       process at most N events\n"
    "      --first-event F      start at entry F (for splitting a file into chunks)\n"
    "      --no-copy-objects    don't copy hAttempts etc. (all but one chunk of a split file)\n"
    "      --count              print the number of events in the input and exit\n"
    "  -h, --help               this message\n";
}

// Interaction channels of the cascade that can be scaled or switched off,
// and the GENIE configuration parameter(s) that do it in each model. The
// cascade multiplies the channel's fate fraction by the factor and
// renormalizes, so 0 switches the channel off. (hN's nucleon elastic /
// inelastic and pion inelastic factors are a local GENIE patch; the others
// are stock GENIE parameters.) The mean-free-path entries scale the
// interaction probability per step rather than a branching fraction.
struct Channel {
  const char *name;
  const char *help;
  vector<string> paramsHN, paramsHA;
  bool isBool;   // INUKE-DoCompoundNucleus: F == 0 -> false, else true
};
const vector<Channel> kChannels = {
  {"nucleon.elastic",    "nucleon elastic scattering",             {"FSI-Nucleon-FracElasScale"}, {}, false},
  {"nucleon.inelastic",  "nucleon inelastic scattering",           {"FSI-Nucleon-FracInelScale"}, {"FSI-Nucleon-FracInelScale"}, false},
  {"nucleon.cex",        "nucleon charge exchange",                {}, {"FSI-Nucleon-FracCExScale"}, false},
  {"nucleon.absorption", "nucleon absorption",                     {}, {"FSI-Nucleon-FracAbsScale"}, false},
  {"nucleon.piprod",     "pion production by nucleons",            {}, {"FSI-Nucleon-FracPiProdScale"}, false},
  {"nucleon.compound",   "compound-nucleus (pre-equilibrium) fate", {"INUKE-DoCompoundNucleus"}, {"INUKE-DoCompoundNucleus"}, true},
  {"pion.elastic",       "pion elastic scattering",                {"INUKE-NucQEFac"}, {}, false},
  {"pion.inelastic",     "pion inelastic scattering",              {"FSI-Pion-FracInelScale"}, {"FSI-Pion-FracInelScale"}, false},
  {"pion.cex",           "pion charge exchange",                   {"INUKE-NucCEXFac"}, {"FSI-Pion-FracCExScale"}, false},
  {"pion.absorption",    "pion absorption",                        {"INUKE-NucAbsFac"}, {"FSI-ChargedPion-FracAbsScale", "FSI-NeutralPion-FracAbsScale"}, false},
  {"pion.piprod",        "pion production by pions",               {}, {"FSI-Pion-FracPiProdScale"}, false},
  {"nucleon.mfp",        "nucleon mean free path (a scale, not a switch)", {"FSI-Nucleon-MFPScale"}, {"FSI-Nucleon-MFPScale"}, false},
  {"pion.mfp",           "pion mean free path (a scale, not a switch)",    {"FSI-ChargedPion-MFPScale", "FSI-NeutralPion-MFPScale"}, {"FSI-ChargedPion-MFPScale", "FSI-NeutralPion-MFPScale"}, false},
};

void ListChannels()
{
  printf("%-20s %-4s %-4s %s\n", "channel", "hN", "hA", "meaning");
  for (const Channel &c : kChannels)
    printf("%-20s %-4s %-4s %s\n", c.name, c.paramsHN.empty() ? "-" : "yes",
           c.paramsHA.empty() ? "-" : "yes", c.help);
  printf("\nUse --channel NAME=FACTOR (0 switches the channel off). A channel marked '-'\n"
         "is not a separate fate in that cascade model.\n");
}

// Write <dir>/{HN,HA}Intranuke2018.xml: copies of the files GENIE would
// otherwise read (the first on the resolved GXMLPATH, so a variant such as
// genie_override_srctail keeps its medium) with the requested parameters set
// in the "Default" parameter set. Returns false with a message on failure.
bool WriteOverride(const string &dir, FSIModel model,
                   const map<string, string> &params, bool quietCascade, string &err)
{
  if (mkdir(dir.c_str(), 0775) != 0 && errno != EEXIST) {
    err = "cannot create " + dir; return false;
  }
  for (const char *base : {"HNIntranuke2018.xml", "HAIntranuke2018.xml"}) {
    const string srcDir = fsi::FindOnGenieXMLPath(base);
    if (srcDir.empty()) { err = string("no ") + base + " on GXMLPATH"; return false; }
    ifstream in(srcDir + "/" + base);
    if (!in) { err = "cannot read " + srcDir + "/" + base; return false; }
    const bool edit = (model == kHN2018) == (string(base) == "HNIntranuke2018.xml");
    ofstream out(dir + "/" + base);
    string line;
    bool inserted = false;
    while (getline(in, line)) {
      if (edit) {
        bool drop = false;
        for (const auto &kv : params)
          if (line.find("name=\"" + kv.first + "\"") != string::npos && line.find("<param") != string::npos)
            drop = true;   // our value replaces it (a commented-out copy is harmless to drop too)
        if (drop) continue;
      }
      out << line << "\n";
      if (edit && !inserted && line.find("<param_set name=\"Default\">") != string::npos) {
        out << "    <!-- set by gcf_fsi (channel option) -->\n";
        for (const auto &kv : params) {
          const bool isBool = kv.second == "true" || kv.second == "false";
          out << "    <param type=\"" << (isBool ? "bool" : "double") << "\" name=\"" << kv.first
              << "\"> " << kv.second << " </param>\n";
        }
        inserted = true;
      }
    }
    if (edit && !inserted) { err = string("no Default parameter set in ") + base; return false; }
    out.close();
    xmlDocPtr doc = xmlParseFile((dir + "/" + base).c_str());
    if (!doc) { err = string("the generated ") + base + " is not well-formed XML"; return false; }
    xmlFreeDoc(doc);
  }
  // With a channel switched off, hadrons left without any allowed fate make
  // INTRANUKE log "Couldn't select a fate" at ERROR level for every such
  // step: expected here, and enough to flood a log. Raise the two cascade
  // streams to FATAL in a copy of the Messenger.xml in use.
  if (quietCascade) {
    const string srcDir = fsi::FindOnGenieXMLPath("Messenger.xml");
    if (srcDir.empty()) { err = "no Messenger.xml on GXMLPATH"; return false; }
    ifstream in(srcDir + "/Messenger.xml");
    if (!in) { err = "cannot read " + srcDir + "/Messenger.xml"; return false; }
    ofstream out(dir + "/Messenger.xml");
    string line;
    while (getline(in, line)) {
      if (line.find("msgstream=\"HNIntranuke2018\"") != string::npos ||
          line.find("msgstream=\"HAIntranuke2018\"") != string::npos) {
        const size_t a = line.find('>'), b = line.find("</priority>");
        if (a != string::npos && b != string::npos && b > a)
          line = line.substr(0, a + 1) + " FATAL " + line.substr(b);
      }
      out << line << "\n";
    }
  }
  return true;
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
  long long maxEvents = 0, firstEvent = 0;
  int defaultA = -1, defaultZ = -1;
  bool independent = false, overwrite = false, copyObjects = true, countOnly = false;
  vector<pair<string, double>> channels;

  enum { kIndependent = 1000, kTarget, kVertexPower, kKeep, kOverwrite, kFirst,
         kNoCopy, kCount, kChannel, kListChannels };
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
    {"first-event", required_argument, nullptr, kFirst},
    {"no-copy-objects", no_argument, nullptr, kNoCopy},
    {"count", no_argument, nullptr, kCount},
    {"channel", required_argument, nullptr, kChannel},
    {"list-channels", no_argument, nullptr, kListChannels},
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
      case kFirst: firstEvent = atoll(optarg); break;
      case kNoCopy: copyObjects = false; break;
      case kCount: countOnly = true; break;
      case kChannel: {
        const char *eq = strchr(optarg, '=');
        char *end = nullptr;
        const double f = eq ? strtod(eq + 1, &end) : -1.;
        if (!eq || end == eq + 1 || *end != '\0' || f < 0.) {
          cerr << "Invalid --channel (expect NAME=FACTOR with FACTOR >= 0, e.g. pion.absorption=0).\n";
          return 1;
        }
        channels.emplace_back(string(optarg, eq - optarg), f);
        break;
      }
      case kListChannels: ListChannels(); return 0;
      default: Usage(); return 1;
    }
  }

  Runner runner;
  bool ok = true;
  auto fail = [&ok](const string &msg) { cerr << "Error: " << msg << "\n"; ok = false; };
  if (input.empty() || (output.empty() && !countOnly)) fail("--input and --output are required");
  if (firstEvent < 0 || maxEvents < 0) fail("--first-event and --max-events must be >= 0");
  if (!input.empty() && input == output) fail("--output must differ from --input");
  if (modelName == "hN") runner.model = kHN2018;
  else if (modelName == "hA") runner.model = kHA2018;
  else fail("--fsi-model must be hN or hA");
  if (vertexPowerName == "auto") runner.vertexPower = fsi::kVertexAuto;
  else if (vertexPowerName == "1" || vertexPowerName == "2" || vertexPowerName == "3")
    runner.vertexPower = atoi(vertexPowerName.c_str());
  else fail("--vertex-power must be auto, 1, 2 or 3");
  map<string, string> overrideParams;
  {
    map<string, double> byName;
    for (const auto &c : channels) {
      const Channel *ch = nullptr;
      for (const Channel &k : kChannels) if (c.first == k.name) ch = &k;
      if (!ch) { fail("unknown channel '" + c.first + "' (see --list-channels)"); continue; }
      const vector<string> &params = (modelName == "hA") ? ch->paramsHA : ch->paramsHN;
      if (params.empty()) {
        fail("channel '" + c.first + "' is not a separate fate in the " + modelName +
             " cascade (see --list-channels)");
        continue;
      }
      byName[c.first] = c.second;
      for (const string &p : params)
        overrideParams[p] = ch->isBool ? (c.second == 0. ? "false" : "true") : to_string(c.second);
    }
    if (modelName == "hA") {
      // GENIE renormalizes these four by 1/sum: all off would divide by zero.
      double sum = 0.;
      for (const char *n : {"nucleon.inelastic", "nucleon.cex", "nucleon.absorption", "nucleon.piprod"})
        sum += byName.count(n) ? byName[n] : 1.;
      if (sum == 0.) fail("hA: at least one of nucleon.inelastic/cex/absorption/piprod must stay on");
    }
  }
  regex keepRe;
  try { keepRe = regex(keep); } catch (const regex_error &) { fail("--keep is not a valid regex"); }
  if (!countOnly && !getenv("GENIE")) fail("GENIE is not set; source setup_env.sh (or build/gcf_env.sh) first");
  if (!ok) { cerr << "(gcf_fsi --help for usage)\n"; return 1; }
  runner.mode = independent ? fsi::kIndependent : fsi::kSharedRemnant;
  if (runner.vertexPower == 2) fsi::SetSRCVertexRhoPower(2);
  // Channel settings go into a private GENIE override directory, first on
  // GXMLPATH, so GENIE reads the modified HN/HAIntranuke2018.xml. Its files
  // are copies of the ones the resolved path would otherwise give, and that
  // path stays behind it, so nothing else changes.
  const string resolvedPath = fsi::ResolveGenieXMLPath();
  string overrideDir;
  if (!overrideParams.empty()) {
    overrideDir = output + ".genie_override";
    bool anyOff = false;
    for (const auto &c : channels) if (c.second == 0.) anyOff = true;
    if (anyOff)
      cout << "Note: with a channel switched off, GENIE's per-step 'Couldn't select a fate' "
              "messages are expected and are suppressed (cascade log streams set to FATAL).\n";
    string err;
    if (!WriteOverride(overrideDir, runner.model, overrideParams, anyOff, err)) {
      cerr << "Error: " << err << "\n"; return 1;
    }
    setenv("GXMLPATH", (overrideDir + ":" + resolvedPath).c_str(), 1);
    fsi::ResolveGenieXMLPath();
  }
  if (!countOnly) fsi::PrintGenieXMLConfig(cout);

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
  if (countOnly) { cout << "Entries: " << nEntries << "\n"; return 0; }
  const long long begin = min(firstEvent, nEntries);
  const long long end = (maxEvents > 0) ? min(nEntries, begin + maxEvents) : nEntries;
  runner.total = end - begin;
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
  if (begin > 0 || end < nEntries) node = node.Range(begin, end);

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
      if (!copyObjects || name == treeName || name == "gcf_fsi" || copied.count(name)) continue;
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
         << "channels:" << (channels.empty() ? " default" : "");
    for (const auto &c : channels) prov << " " << c.first << "=" << c.second;
    prov << "\n"
         << "events: " << runner.done << ", cascade ran: " << runner.nRan
         << ", vertex sampled: " << runner.nVertexSampled << "\n";
    out.cd();
    TNamed("gcf_fsi", prov.str().c_str()).Write("gcf_fsi", TObject::kOverwrite);
  }

  if (!overrideDir.empty()) {   // the settings are in the provenance record
    for (const char *base : {"HNIntranuke2018.xml", "HAIntranuke2018.xml", "Messenger.xml"})
      remove((overrideDir + "/" + base).c_str());
    rmdir(overrideDir.c_str());
  }
  cout << "FSI (" << modelName << ", " << (independent ? "independent" : "shared")
       << ") applied to " << runner.done << " events in " << timer.RealTime() << " s: "
       << "cascade ran for " << runner.nRan << ", skipped for "
       << runner.done - runner.nRan << " (residual A < 3 or nothing to transport); "
       << runner.nVertexSampled << " vertices sampled.\nWrote " << output << "\n";
  return 0;
}
