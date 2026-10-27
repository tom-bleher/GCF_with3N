# genQE_FSI: GCF Quasi-Elastic Generator with GENIE Final State Interactions

## Overview

This project generates quasi-elastic (QE) electron-nucleus scattering events using the **Generalized Contact Formalism (GCF)** for Short-Range Correlations (SRC), with nuclear **Final State Interactions (FSI)** modeled by GENIE's intranuclear cascade.

The physics pipeline is:

```
e + A  -->  e' + N_lead + N_recoil + (A-2) residual
                  |           |
                  +--- GENIE INTRANUKE FSI ---+
                  |           |
              N'_lead     N'_recoil  (+secondaries: pions, knocked-out nucleons, ...)
```

The GCF generator produces an SRC nucleon pair (lead + recoil) with correlated momenta from a contact-formalism wavefunction. Each nucleon then propagates through the (A-2) residual nucleus via GENIE's INTRANUKE cascade, which can elastically scatter, charge-exchange (p <-> n), produce pions, or absorb the nucleon. Pauli blocking rejects events where the final nucleon momentum falls below the Fermi surface.

### FSI Models

Two GENIE cascade models are available, selectable at runtime:

| Model | Class | Description | Speed |
|---|---|---|---|
| **hN** (default) | `HNIntranuke2018` | Full intranuclear cascade. Steps the hadron through the nucleus in 0.05 fm increments, resolving individual hadron-nucleon collisions. Includes medium corrections (Salcedo-Oset for pions, [doi:10.1016/0375-9474(88)90310-7](https://doi.org/10.1016/0375-9474(88)90310-7); Pandharipande-Pieper for NN, [doi:10.1103/PhysRevC.45.791](https://doi.org/10.1103/PhysRevC.45.791)). | ~1x |
| **hA** | `HAIntranuke2018` | Effective cascade. Uses total hadron-nucleus cross sections to decide the outcome in one step. Faster but less detailed. | ~5-10x |

For SRC studies, **hN is recommended** because it accurately models nucleon rescattering with proper angular distributions and secondary particle production. Use hA for fast exploratory runs or systematics scans.

---

## Prerequisites

| Dependency | Version | Install |
|---|---|---|
| ROOT | 6.x | `brew install root` or from source |
| GSL | 2.x | `brew install gsl` |
| log4cpp | any | `brew install log4cpp` |
| libxml2 | any | `brew install libxml2` |
| CMake | >= 3.6 | `brew install cmake` |
| Python 3 | 3.x | For plotting scripts (matplotlib, numpy) |

GENIE R-3_06_02 is built locally inside this directory (see below). PYTHIA6 and LHAPDF are **not** required.

---

## Building

### Step 1: Build GENIE (one-time setup)

Clone and build GENIE R-3_06_02 inside this project directory:

```bash
cd src/programs/genQE_FSI

git clone --depth 1 --branch R-3_06_02 \
  https://github.com/GENIE-MC/Generator.git Generator-R-3_06_02

cd Generator-R-3_06_02
export GENIE=$(pwd)

./configure \
  --prefix=$(pwd)/install \
  --disable-pythia6 \
  --disable-lhapdf5 \
  --disable-flux-drivers \
  --disable-geom-drivers \
  --with-compiler=clang \
  --with-libxml2-inc=/opt/homebrew/opt/libxml2/include/libxml2 \
  --with-libxml2-lib=/opt/homebrew/opt/libxml2/lib \
  --with-log4cpp-inc=/opt/homebrew/opt/log4cpp/include \
  --with-log4cpp-lib=/opt/homebrew/opt/log4cpp/lib

make -j$(sysctl -n hw.ncpu)
```

> **Note:** If the HEDIS module fails to compile (LHAPDF dependency), remove its lines from the top-level Makefile or build only the needed targets:
> ```bash
> make framework physics-utilities physics-nuclear-environment \
>      physics-hadronic-simulations physics-neutrino-scattering-modes
> ```

### Step 2: Build the project

```bash
cd src/programs/genQE_FSI
mkdir -p build && cd build
cmake ..
make -j$(sysctl -n hw.ncpu)
```

This produces three executables in `build/`:

| Binary | Purpose |
|---|---|
| `genQE` | QE generator **without** FSI (baseline) |
| `genQE_FSI` | QE generator **with** FSI, outputs ROOT tree |
| `SRC_analysis_2N` | Analysis program: generates events and fills histograms |

CMake automatically embeds the GENIE library path via `RPATH`, so no `DYLD_LIBRARY_PATH` is needed at runtime.

To build without FSI (no GENIE dependency):

```bash
cmake .. -DUSE_FSI=OFF
make -j$(sysctl -n hw.ncpu)
```

---

## Running

### Quick start with `run_SRC.sh`

The wrapper script sets all GENIE environment variables and launches `SRC_analysis_2N`:

```bash
cd src/programs/genQE_FSI

# 100k events (default, hN model)
./run_SRC.sh

# Custom event count
./run_SRC.sh 500000

# Disable FSI (pass 0 as second argument)
./run_SRC.sh 100000 0

# Use hA model instead of hN
./run_SRC.sh 100000 1 hA
```

### Running `SRC_analysis_2N` directly

```bash
export GENIE="$(pwd)/Generator-R-3_06_02"
export GXMLPATH="$(pwd)/config/genie_override:$GENIE/config/G18_10a:$GENIE/config"
export GMSGCONF="$(pwd)/config/quiet_messenger.xml"

./build/SRC_analysis_2N <nEvents> [doFSI] [model]
```

Arguments:
- `nEvents`: Number of accepted events to generate (default: 1000000)
- `doFSI`: 1 = FSI enabled (default), 0 = FSI disabled
- `model`: `hN` (default) or `hA` — selects the GENIE cascade model

The program runs on **carbon-12** with the **on-shell** (Rosenbluth) eN cross section, both hardcoded. Later positional arguments set σ_CM, the sampling mode (2N pair, SRC-only or single-nucleon mean field), the beam energy (default 5.01 GeV) and more; the full list is in the header of `analysis/SRC_analysis_2N.cpp`.

The Q² sampling range defaults to 1–10 GeV². Set `GCF_QSQ_RANGE=MIN:MAX` (GeV²) to change it, e.g. `GCF_QSQ_RANGE=3:25` for EIC-energy runs. The range is printed at startup and written to the `.meta.txt` file.

### Running `genQE_FSI` (ROOT tree output)

For custom targets, beam energies, and ROOT tree output:

```bash
export GENIE="$(pwd)/Generator-R-3_06_02"
export GXMLPATH="$(pwd)/config/genie_override:$GENIE/config/G18_10a:$GENIE/config"
export GMSGCONF="$(pwd)/config/quiet_messenger.xml"

./build/genQE_FSI <Z> <N> <Ebeam_GeV> <output.root> <nEvents> [flags]
```

Example (carbon-12, 4.5 GeV, 500k events):

```bash
./build/genQE_FSI 6 6 4.5 output_C12.root 500000
```

**Flags:**
| Flag | Description |
|---|---|
| `-v` | Verbose output (print progress every 100k events) |
| `-u <type>` | NN interaction potential (default: `AV18`, auto `AV18_deut` for deuterium) |
| `-c <model>` | eN cross section model: `cc1`, `cc2`, or `onshell` |
| `-s <sigma>` | Center-of-mass width sigma_CM [GeV/c] |
| `-E <Estar>` | Separation energy E* [GeV] |
| `-M` | Use randomized E* (Barack's values) |
| `-k <kCut>` | Relative momentum hard cutoff [GeV/c] (default: 0.25) |
| `-O` | Enable peaking radiation |
| `-C` | Enable Coulomb correction |
| `-r` | Randomize nuclear properties |
| `-l` | Use lightcone cross section |

**ROOT tree branches** (tree name: `genT`):
| Branch | Type | Description |
|---|---|---|
| `lead_type` | `Int_t` | PDG-like code for lead nucleon (2212=p, 2112=n) |
| `rec_type` | `Int_t` | PDG-like code for recoil nucleon |
| `pe[3]` | `Double_t` | Scattered electron 3-momentum [GeV/c] |
| `pLead[3]` | `Double_t` | Lead nucleon 3-momentum after FSI [GeV/c] |
| `pRec[3]` | `Double_t` | Recoil nucleon 3-momentum after FSI [GeV/c] |
| `pAm2[3]` | `Double_t` | (A-2) residual momentum [GeV/c] |
| `pRel[3]` | `Double_t` | Pair relative momentum [GeV/c] |
| `q[3]` | `Double_t` | Virtual photon 3-momentum [GeV/c] |
| `weight` | `Double_t` | Event weight (0 if absorbed/Pauli-blocked) |

### EIC (collider-frame) output

The generators work in the target rest frame with the electron along +z. An EIC collision of an electron of energy E_e with an ion of momentum A·p_N is the same event in the rest frame with an electron of energy E_IRF = γ(1+β)E_e, where γβ = A·p_N/M_A. Run the generator at E_IRF and convert its output with `tools/gcf2hepmc.py`. The converter rotates by π about y (electron along −z), boosts with the ion velocity, and writes head-on HepMC3 with the electron and the ion as status-4 beams. Feed the result to the afterburner (crossing angle, beam effects) and then npsim.

```bash
# rest-frame electron energy: 440.33 GeV for e+12C at 5x41, 2362.49 GeV at 10x110
E=$(../../../tools/gcf2hepmc.py --beams 5x41 --ion C12 --print-ebeam-irf)
printf 'QSq 3 25\n' > ps_eic.txt
./build/genQE 6 6 $E gen.root 200000 -P ps_eic.txt -M          # PRC 105, 034001 settings (doi:10.1103/PhysRevC.105.034001)
../../../tools/gcf2hepmc.py gen.root gen.hepmc3 --beams 5x41 --ion C12 --paper-cuts
# inside eic-shell:
abconv -p ip6_eAu_41x5 -f hepmc3 -o ab gen.hepmc3
npsim --compactFile geometry.xml --inputFiles ab.hepmc --outputFile sim.edm4hep.root -N -1
```

- **Production driver.** `tools/run_eic_sim.sh` runs the whole chain: genQE on the host, then abconv, chunked npsim and EICrecon inside eic-shell. Each run directory gets a `run.json` recording versions, geometry, seeds, material map and job status. For example:
  ```bash
  EPIC_INSTALL=~/eic/dev/epic_craterlake_5x41/install EICRECON_INSTALL=~/eic/dev/EICrecon/install \
    tools/run_eic_sim.sh --beams 5x41 --ion C12 --attempts 1500000 --seed 1 --out runs/c12_5x41
  ```
  - Default geometry: `epic_craterlake_5x41` with `fields/beamline_5x41_He4.xml` (Z/A = ½ rigidity), and `epic_craterlake_10x110_H2` for 10x110.
  - abconv always draws the same random sequence (it has no seed option), so the script runs it once over the whole sample and only then splits the events into jobs.
  - **Q² slices.** With flat Q² sampling, Q² = 3–5 GeV² carries about 90% of the cross section but only about 9% of the events. Production therefore runs disjoint slices (`--q2 3:5`, `--q2 5:12`, `--q2 12:25`), each a separate run with its own seed. Each slice is normalized by its own attempts (in `run.json`), and the slices add: 0.0312 + 0.0039 + 0.0001 nb equals the single 3–25 GeV² run, 0.035 nb (¹²C, 5x41, paper cuts).
  - **Other generators.** `--generator genQE_FSI` runs GENIE FSI. Add `--gen-args "-n"` to switch FSI off (a PWIA baseline from the same generator; it reproduces genQE to 2%), or `--gen-args "-n -F 0.25" --no-cuts` for single-nucleon mean-field knockout (PRC 105, 034001 Fig. 7, [doi:10.1103/PhysRevC.105.034001](https://doi.org/10.1103/PhysRevC.105.034001)). E* sampling (`-M`) is added only for SRC pairs in A > 2 nuclei: for the deuteron it would give the empty remnant a mass.
  - **Beyond quasielastic 2N.** `--generator genQE_3N_FSI` (3N SRC; `--attempts` counts accepted events, about 8% pass the paper cuts; the rest-frame electron window scales as 1/E_IRF) and `--generator genDIS --no-cuts` (tagged DIS off an SRC pair with CT18NN PDFs; the electron, recoil and A-2 remnant are simulated, the unhadronized struck system is not). The analysis makes the 3N recoil-tagging figure (0, 1 or 2 recoils detected, by recoil pair) and the tagged-DIS figures (recoil acceptance vs the spectator light-cone fraction alpha_s and momentum; tagged yield over x' and Q2).
  - **Reconstructed momenta.** The analysis stores each nucleon's reconstructed momentum (B0/central track, or reconstructed ZDC neutron) boosted back to the ion rest frame with the nominal crossing angle and boost. Fig. 2 gives dp/p vs p; Fig. 6's blue points are the spectrum in reconstructed p_IRF, the analog of the paper's smeared points. RP/OMD have no reconstruction for ion beams.
  - **Evaporation.** `--deexcite` breaks up the excited A−2 remnant (E* from genQE `-M`) before the afterburner, with one of two back ends:
    - `--deexcite` or `--deexcite geant4`: Geant4's nuclear de-excitation (evaporation, Fermi break-up, photon emission), `tools/deexcite` (build it inside eic-shell: `cmake -S tools/deexcite -B tools/deexcite/build && cmake --build tools/deexcite/build`). It breaks the remnant up at rest with mass M_gs + E* and boosts the products back.
    - `--deexcite beagle`: BeAGLE with FLUKA, the evaporation of PRC 105, 034001 itself, `tools/beagle_deexcite.py`. Each event goes to BeAGLE's `MODEL GCF-FT` mode at nuclear-reaction level 2, which hands the remnant to FLUKA and rescatters nothing. It needs BeAGLE built against CERN FLUKA 4 (`BEAGLESYS`, `FLUPRO`, `LHAPATH`; see BeAGLE's `build_local.sh`). FLUKA conserves the remnant's energy to the keV level BeAGLE writes and its momentum to about 10⁻⁴ of itself (up to ~0.2 MeV/c).

    Both check A, Z and four-momentum conservation, keep the remnant as a status 3 entry and append the products after all other particles, so the analysis finds e', lead and recoil where it always does. `gen_deex.hepmc3.json` records the products and checks. With the seed of an earlier genQE run and fewer `--attempts`, the generator writes exactly that run's first events, so an evaporation run can be compared event by event with a run without it. The analysis then adds `figEvap_angles_momenta.pdf` (Fig. 3 with evaporation nucleons, their rest-frame momenta) and an `evaporation` block in `summary.json` (multiplicities, far-forward detection, overlap with tagged recoil neutrons in the ZDC).
  - **Contact uncertainty.** `tools/eic_src_analysis.py band` varies the ¹²C contacts (AV18, k-space, with their uncertainties from `gcfSRC.cc`, independently as `genQE -r` does) and reweights every simulated event exactly by S'(k_rel)/S(k_rel), using the generator's `pRel`. The weight is linear in the contacts, so the band costs four weighted sums per observable. It writes the Q² yield, the accepted recoil spectra and the relative uncertainty with the 68% band, and percentiles of σ, tagged yields, acceptances and pair shares in `summary.json`.
  - **Compression.** Jobs write to local scratch and land a ZSTD-5 copy (`--compress 505`, the default) in the run directory. That is about a third smaller than npsim's and EICrecon's ZLIB-1 output, reads faster, and holds identical contents.
  - **Analysis.** `tools/eic_src_analysis.py` handles FSI output. The leading or recoil nucleon is the highest-momentum nucleon descended from the generated one, the detector follows its final type, and absorbed nucleons count as undetected. Mean-field runs give Fig. 7 (p_miss < 0.3 GeV/c, E_miss < 80 MeV, Q² > 3, electron and nucleon detected).
  - **Resume.** `--stage resume` finishes an interrupted run. It reruns npsim and EICrecon only for the chunks without a successful exit code.
  - **Long runs.** Start them detached so they survive the terminal, e.g. `setsid nohup tools/run_eic_sim.sh ... > run.out 2>&1 &`.
- **Seeds.** `genQE -S <seed>` makes a run reproducible. The seed is printed and stored in the output (`seed`), and the converter copies it into its JSON. With no `-S`, a random seed is drawn and printed.
- **Which genQE.** Use this directory's `build/genQE`: its σ_CM default for ¹²C is 150 MeV/c (as in PRC 105, 034001), whereas the suite's `src/programs/genQE` defaults to 55 MeV/c. `-M` samples E* as in that paper (17.3 ± 9.6 MeV).
- **What gets written.**
  - `genQE` events: the scattered electron, the leading and recoil nucleons, and the A−2 remnant nucleus (mass from the `mAm2` branch, including E*). An unbound remnant (nn, pp) is skipped.
  - `SRC_analysis_2N`, `genQE_FSI` and `genQE_3N_FSI` events: the electron plus the final-state particle list. GENIE drops the residual nucleus, so these have no remnant.
- **Bookkeeping.** The HepMC event number is the source tree entry; the event weight is the generator weight (nb per attempt). `<output>.json` records the beams, cuts, counts, attempts and the cross section of the written sample.
- **Beam check.** The converter refuses a file generated at the wrong E_IRF. Trees without `q` (`genQE_3N_FSI`) need `--ebeam-irf`.
- **genQE_3N_FSI.** The default electron-angle window (`-t 5:40` degrees) is for fixed-target energies. Use `-t 0.1:1.5 -q 3:25` at 5x41 and `-t 0.02:0.3 -q 3:25` at 10x110.
- **Remnant hits.** The remnant nucleus is tracked by Geant4. A ¹⁰Be or ¹⁰C remnant (Z/A ≠ ½) reaches the Roman pots, so attribute far-forward hits to particles rather than counting them per event.

---

## Output

### `SRC_analysis_2N` output

Results are written to `analysis_output_2N/txt_files/`. Key files:

| File | Contents |
|---|---|
| `hist_theta12_1D.txt` | Angle between lead and recoil nucleons (100 bins, 0-180 deg) |
| `hist_xB_Q2.txt` | 2D Bjorken-x vs Q^2 distribution (500x500 bins) |
| `hist_theta12_theta23_3body.txt` | 2D theta heatmap with hypothetical 3rd nucleon |
| `hist_<varname>_1D.txt` | 1D histograms for all kinematic variables (45 bins each) |
| `hist_<varname>_region<N>_1D.txt` | Same histograms split by xB region |
| `hist_xB_Q2_region<N>.txt` | 2D xB-Q2 per region |
| `hist_theta12_theta23_3body_region<N>.txt` | 2D theta heatmap per region |
| `fsi_event_stats.txt` | FSI channel summary (charge exchange, elastic, pion rates) |

**xB regions:**
| Region | Range | Description |
|---|---|---|
| 0 | xB < 1.0 | Quasi-elastic peak |
| 1 | 1.0 < xB < 1.5 | 2N SRC onset |
| 2 | 1.5 < xB < 2.0 | 2N SRC region |
| 3 | 2.0 < xB < 2.5 | High xB |
| 4 | xB > 2.5 | Very high xB |

**Kinematic variables saved** (1D histograms):
Q2, xB, outgoing electron angle/momentum, incoming lead angle/pmiss, recoil angle/momentum, outgoing lead angle/momentum, theta12, pair CM momentum, pair relative momentum, hypothetical 3rd nucleon momentum/angles.

### Console output

The program prints a summary after completion:

```
--- FSI / event statistics ---
  FSI mode:             ENABLED (GENIE intranuke transport)
  Total attempts:       ...
  Accepted (weight>0):  ...
  Gen. efficiency:      ~57%

--- FSI channel summary (accepted events) ---
  Events with any CX:    ... (~4%)
  Events with elastic:   ... (~1%)
  Events with any pion:  ... (~4%)
  Total pions:           ... [pi+=.., pi-=.., pi0=..]
```

---

## Plotting

Two Python scripts are provided in `plotting/`:

### 1D histograms

```bash
cd src/programs/genQE_FSI
python3 plotting/plot_SRC_analysis_2N.py
```

Reads all `hist_*_1D.txt` files from `analysis_output_2N/txt_files/` and produces PNG plots in `analysis_output_2N/png_files/`. Automatically loads FSI statistics from `fsi_event_stats.txt` if present.

### 2D theta heatmaps

```bash
python3 plotting/plot_theta_heatmap_2N.py
```

Plots 2D theta12-theta23 heatmaps (lead-recoil angle vs recoil-hypothetical-3rd-nucleon angle), both for all events and per-region.

---

## Configuration

### GENIE cascade tuning

XML override files in `config/genie_override/` control GENIE cascade parameters:

- **`HNIntranuke2018.xml`** — hN model parameters
- **`HAIntranuke2018.xml`** — hA model parameters
- **`Messenger.xml`** — Log verbosity (ERROR-only for production)

The programs set `GXMLPATH` themselves (so leaving it unset is fine), in this order: your own `GXMLPATH` entries, then the project override directory, then `$GENIE/config/G18_10a`, then `$GENIE/config`. GENIE reads each file from the first directory that has it, so a directory you put on `GXMLPATH` always wins over the project override. The project override is `config/genie_override`, or the directory in `$GENIE_OVERRIDE_DIR` if set (e.g. `config/genie_override_srctail`); when the variable selects a different directory, a `GXMLPATH` entry naming the default `config/genie_override` (as `build/gcf_env.sh` and the `gcf` front end export) is dropped, so the selection always takes effect. A selected directory without `HNIntranuke2018.xml` gives a warning at startup (`gcf_fsi` stops). The resolved path, and which directory and `NuclearModel` each Intranuke file comes from, are printed at startup.

Key tunable parameters (both models share these via the `Intranuke2018` base class):

| Parameter | Default | Description |
|---|---|---|
| `INUKE-HadStep` | 0.05 fm | Step size for cascade propagation |
| `INUKE-NucRemovalE` | 0.00 GeV | Binding energy subtracted from cascade nucleons |
| `INUKE-FermiMomentum` | 0.250 GeV/c | Fermi momentum for nuclear medium |
| `INUKE-NucAbsFac` | 1.0 | Absorption cross section scale factor |
| `INUKE-NucQEFac` | 1.0 | Quasi-elastic cross section scale factor |
| `INUKE-NucCEXFac` | 1.0 | Charge exchange cross section scale factor |
| `INUKE-XsecNNCorr` | true | Pandharipande-Pieper medium corrections for NN ([doi:10.1103/PhysRevC.45.791](https://doi.org/10.1103/PhysRevC.45.791)) |
| `INUKE-DoFermi` | true | Enable Fermi motion of target nucleons |
| `INUKE-DoCompoundNucleus` | true | Enable compound nucleus formation |

hN-specific parameters:

| Parameter | Default | Description |
|---|---|---|
| `HNINUKE-UseOset` | true | Salcedo-Oset medium corrections for pion-nucleon ([doi:10.1016/0375-9474(88)90310-7](https://doi.org/10.1016/0375-9474(88)90310-7)) |
| `HNINUKE-AltOset` | false | Alternative Oset parameterization |

### Fermi momentum (Pauli blocking)

The Fermi momentum for Pauli blocking is set in the analysis code:

```cpp
myGen->SetFSITuning(220.);  // MeV/c
```

Typical values: 220 MeV/c (carbon), 250 MeV/c (iron), 270 MeV/c (lead).

### Log suppression

`config/quiet_messenger.xml` silences all GENIE log output below ERROR level. Without this, GENIE produces ~1 GB of log output per 10k events. The `GMSGCONF` environment variable must point to this file.

---

## Using QEGeneratorFSI in your own code

### Basic usage

```cpp
#include "QEGeneratorFSI.hh"

// Initialize
gcfNucleus* nucleus = new gcfNucleus(6, 6, (char*)"AV18");  // Carbon-12
eNCrossSection* cs = new eNCrossSection(cc1, kelly);
TRandom3* rng = new TRandom3(0);
QEGeneratorFSI* gen = new QEGeneratorFSI(4.5, nucleus, cs, rng);

// Configure FSI
gen->EnableFSI(true);                // on by default
gen->SetFSIModel(kHN2018);          // or kHA2018 for fast mode
gen->SetFSITuning(220.);            // Fermi momentum in MeV/c

// Event loop
for (int i = 0; i < nEvents; i++) {
    double weight;
    int lead_type, rec_type;
    TLorentzVector vk, vLead, vRec, vAm2;
    TVector3 vRel;
    TLorentzVector vq;
    double Estar;

    gen->generate_event(weight, lead_type, rec_type,
                        vk, vLead, vRec, vAm2, vRel, vq, Estar);

    if (weight <= 0.) continue;  // absorbed or Pauli-blocked

    // Use vLead, vRec (post-FSI momenta)
}
```

### Accessing FSI details

After each event, you can inspect what happened during FSI:

```cpp
// Per-event FSI summary flags
const auto& stats = gen->GetLastFSIEventStats();
if (stats.leadChargeExchange)  { /* lead nucleon changed identity */ }
if (stats.recoilElasticLike)   { /* recoil scattered elastically */ }
if (stats.nPionsTotal > 0)     { /* pions were produced */ }

// Full list of FSI secondaries (pions, knocked-out nucleons, etc.)
const auto& secondaries = gen->GetLastFSISecondaries();
for (const auto& sec : secondaries) {
    // sec.pdg       — PDG code (211=pi+, -211=pi-, 111=pi0, 2212=p, 2112=n)
    // sec.p4        — 4-momentum (TLorentzVector)
    // sec.parentRole — 0=from lead transport, 1=from recoil transport
}
```

### Switching FSI models

```cpp
// Use hN (default, more accurate)
gen->SetFSIModel(kHN2018);

// Use hA (faster, for systematics)
gen->SetFSIModel(kHA2018);
```

### Comparing FSI on vs off

```cpp
// Run without FSI
gen->EnableFSI(false);
for (int i = 0; i < nEvents; i++) {
    gen->generate_event(weight, lead_type, rec_type, vk, vLead, vRec, vAm2);
    // Fill "no FSI" histograms
}

// Run with FSI
gen->EnableFSI(true);
for (int i = 0; i < nEvents; i++) {
    gen->generate_event(weight, lead_type, rec_type, vk, vLead, vRec, vAm2);
    if (weight <= 0.) continue;
    // Fill "with FSI" histograms
}
```

---

## Verification

Run the integration test suite:

```bash
cd src/programs/genQE_FSI
bash test_fsi_integration.sh
```

This runs 24 checks covering:
1. GENIE library architecture and linkage
2. Build (cmake + make + binary existence)
3. RPATH embedding (no DYLD_LIBRARY_PATH needed)
4. Runtime (1000 events without crash)
5. Physics validation (charge exchange ~2-8%, pion production ~1-10%, efficiency ~40-80%, all pion species present)
6. Output file generation
7. Legacy file removal (no INukeNNData or SimpleFSIHandler)
8. Clean preprocessor guards (no `USE_GENIE_FSI` references)

---

## Project Structure

```
genQE_FSI/
  CMakeLists.txt              # Build system (cmake)
  run_SRC.sh                  # Run wrapper with GENIE env setup
  test_fsi_integration.sh     # Integration test suite (24 checks)
  README_FSI.md               # This file

  QEGeneratorFSI.hh           # FSI generator header (FSIModel enum, API)
  QEGeneratorFSI.cc           # FSI implementation (GENIE cascade integration)
  QEGenerator.hh / .cc        # Base QE generator (no FSI)
  genQE.cc                    # Main for genQE binary (no FSI)
  genQE_FSI.cc                # Main for genQE_FSI binary (ROOT tree output)

  analysis/
    SRC_analysis_2N.cpp        # Histogram-based analysis (deuterium, 6 GeV)

  plotting/
    plot_SRC_analysis_2N.py    # 1D histogram plotter
    plot_theta_heatmap_2N.py   # 2D theta heatmap plotter

  config/
    quiet_messenger.xml        # GENIE log suppression
    genie_override/
      HNIntranuke2018.xml      # hN cascade tuning parameters
      HAIntranuke2018.xml      # hA cascade tuning parameters
      Messenger.xml            # Per-stream GENIE log levels

  Generator-R-3_06_02/         # GENIE source + libs (not tracked in git)
  build/                       # CMake build directory
  analysis_output_2N/          # Output from SRC_analysis_2N
    txt_files/                 # Histogram text files
    png_files/                 # Plot images
```

---

## Light Nuclei Handling

For targets with A <= 4 (e.g., deuterium, He-3, He-4), the residual nucleus after removing the SRC pair has A_res = A - 2:

- **A_res <= 0** (deuterium, A=2): No FSI applied. The SRC pair exits without rescattering.
- **0 < A_res < 7** (He-3, He-4): GENIE FSI runs on the residual, but a warning is printed. GENIE's cascade is validated for A > 6 (heavier than lithium); results for very light residuals should be treated with caution.
- **A_res >= 7** (carbon and heavier): Normal FSI, fully validated.

---

## References

- GENIE Physics and User Manual (v3.06): INTRANUKE models (hA, hN)
- S. Dytman, "FSI models in GENIE" (2025, CERN Indico)
- Papadopoulou et al. (e4nu), Phys. Rev. D 103, 113003 (2021), [doi:10.1103/PhysRevD.103.113003](https://doi.org/10.1103/PhysRevD.103.113003) — GENIE FSI validation with electron scattering
- Harewood & Gran, arXiv:1906.10576 (2019), [doi:10.48550/arXiv.1906.10576](https://doi.org/10.48550/arXiv.1906.10576) — GENIE hA elastic scattering corrections
- Nikolakopoulos et al., Phys. Rev. C 105, 054603 (2022), arXiv:2202.01689, [doi:10.1103/PhysRevC.105.054603](https://doi.org/10.1103/PhysRevC.105.054603) — Cascade model validation
