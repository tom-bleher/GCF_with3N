#!/usr/bin/env bash
# Run genQE quasielastic SRC events through the ePIC chain:
#   genQE (rest frame, on the host) -> tools/gcf2hepmc.py (collider frame)
#   -> abconv (crossing angle, beam effects) -> npsim -> eicrecon   (inside eic-shell)
#
# Usage:
#   tools/run_eic_sim.sh --beams 5x41 --ion C12 --attempts 1500000 --out DIR [options]
#
# Options:
#   --beams 5x41|10x110   electron GeV x ion GeV/c per nucleon (required)
#   --ion C12|d|He4       target; must have Z/A = 1/2, the rigidity of the light-ion lattices (required)
#   --attempts N          genQE attempts; about 13% pass the paper cuts (required)
#   --out DIR             run directory, created; must not exist unless --stage is given (required)
#   --q2 MIN:MAX          Q2 sampling range in GeV^2 (default 3:25). Production can be split into
#                         disjoint Q2 slices, each its own run with its own seed and attempts; a
#                         slice's events normalize by its own attempts, and slices add.
#   --seed S              genQE seed (default 1); npsim chunk i uses S*100000 + i + 1
#   --max-events N        convert at most N events (smoke tests)
#   --chunk-size N        events per npsim/eicrecon job (default 500)
#   --jobs N              parallel npsim jobs (default 64; up to ~2.5 GB each)
#   --reco-jobs N         parallel eicrecon jobs (default 40; up to ~3.1 GB each with one thread). Queue
#                         runs one after another: two at these defaults need more memory than the machine has
#   --reco-threads N      eicrecon threads per job (default 1)
#   --generator G         genQE (default; SRC pairs, no FSI) or genQE_FSI (GENIE FSI; with
#                         --gen-args "-n" no FSI, with "-F 0.25" single-nucleon mean field),
#                         genDIS (DIS off an SRC pair, CT18NN PDFs; use --no-cuts), or genQE_3N_FSI
#                         (3N SRC; --attempts is then the number of accepted events; "-n" for no FSI)
#   --gen-args "ARGS"     extra generator flags. E* sampling (-M) is added for A > 4 SRC runs (4He runs
#                         at E* = 0; Pybus et al. 2020, doi:10.1016/j.physletb.2020.135429, varied the
#                         4He A-2 E* over 0-10 MeV as a systematic)
#   --deexcite [geant4|beagle]
#                         de-excite the A-2 remnant (evaporation nucleons, light fragments, photons) before
#                         the afterburner (genQE and genDIS, whose output carries the remnant):
#                         geant4 (default) uses tools/deexcite, built inside eic-shell with
#                         cmake -S tools/deexcite -B tools/deexcite/build && cmake --build tools/deexcite/build;
#                         beagle uses tools/beagle_deexcite.py, BeAGLE with FLUKA as PRC 105, 034001 did
#                         (needs BeAGLE built against FLUKA 4; see BEAGLESYS, FLUPRO, LHAPATH below)
#   --no-cuts             convert all events instead of Q2 > 3, xB > 1.2
#   --compress N          ROOT compression of sim/reco output (default 505 = ZSTD level 5, about a
#                         third smaller than npsim/EICrecon's ZLIB 1; 0 keeps their output as is)
#   --stage all|gen|sim|reco|resume|record
#                         run one stage only (default all). resume finishes an interrupted run: it
#                         keeps the split chunks and runs npsim and eicrecon only for the chunks
#                         that have not succeeded yet. record rewrites run.json.
#
# Environment:
#   EPIC_INSTALL      epic install prefix (required), e.g. ~/eic/dev/epic_craterlake_5x41/install
#   EICRECON_INSTALL  EICrecon install prefix (required for reco)
#   EIC_SHELL         eic-shell launcher (default ~/eic/eic-shell)
#   BEAGLESYS, FLUPRO, LHAPATH
#                     BeAGLE checkout, FLUKA 4 installation and LHAPDF 5 PDF sets for --deexcite beagle
#                     (defaults ~/src/BeAGLE, ~/src/beagle-deps/fluka, ~/src/beagle-deps/lhapdf5/share/lhapdf/PDFsets)
#   FIELDMAPS         directory linked as ./fieldmaps in every job (default ~/eic/fieldmaps)
#   EPIC_CONFIG       detector config (default epic_craterlake_5x41 for 5x41, epic_craterlake_10x110_H2 for 10x110)
#   FIELD             beamline field file swapped into the config (default fields/beamline_5x41_He4.xml for
#                     5x41; none for 10x110, whose H2 config already has the Z/A = 1/2 field)
#   AB_PRESET         abconv preset (default ip6_hidiv_41x5 / ip6_eAu_110x10). abconv has no light-ion
#                     preset, so C12, d and He4 borrow one, an approximation: 25 mrad crossing angle
#                     plus 0.1 mrad vertical, that preset's ion divergence and bunch length, no
#                     momentum spread. At 5x41 the stand-in is the 41 GeV proton beam (7.5 cm bunch,
#                     vertex z spread about 38 mm rms): ip6_eAu_41x5's 11.6 cm is the Au (and Pb)
#                     value at 41 GeV/u and the longest of any preset, where p, d, 3He, Cu, Ru and Au
#                     at 100-166 GeV/u all use 6-7.5 cm; with it, 17% of 5x41 Q2 12-25 events lose the
#                     electron track (vertex z < -100 mm misses the vertex and Si barrel layers, #190).
#                     PRC 105, 034001 (doi:10.1103/PhysRevC.105.034001) applied only the crossing
#                     angle and parametrized smearing. Names are case sensitive (abconv's help
#                     lists ip6_eau_41x5, which it rejects); ab.hepmc's beam metadata (ion_beam_pdg 11,
#                     ion_beam_energy 18) is wrong, and run.json records the preset used.
#
# The run directory gets gen.root, gen.hepmc3 (+ .json), ab.hepmc, geometry.xml, chunk_*.hepmc,
# sim_*.edm4hep.root, reco_*.edm4eic.root, per-job logs and exit codes, and run.json with the
# software versions, geometry, seeds and counts. abconv runs once over the whole sample (its
# random sequence is fixed, so per-chunk runs would repeat the same beam smearing).
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
GCF=$(cd "$HERE/.." && pwd)

BEAMS="" ION="" ATTEMPTS="" OUT="" Q2="3:25" SEED=1 MAXEV="" CHUNK=500 JOBS=64 RJOBS=40 RTHREADS=1 CUTS="--paper-cuts" STAGE=all GEN=genQE GEN_ARGS="" COMPRESS=505 DEEX=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --beams) BEAMS=$2; shift 2 ;;
    --ion) ION=$2; shift 2 ;;
    --attempts) ATTEMPTS=$2; shift 2 ;;
    --out) OUT=$2; shift 2 ;;
    --q2) Q2=$2; shift 2 ;;
    --seed) SEED=$2; shift 2 ;;
    --max-events) MAXEV=$2; shift 2 ;;
    --chunk-size) CHUNK=$2; shift 2 ;;
    --jobs) JOBS=$2; shift 2 ;;
    --reco-jobs) RJOBS=$2; shift 2 ;;
    --reco-threads) RTHREADS=$2; shift 2 ;;
    --no-cuts) CUTS=""; shift ;;
    --generator) GEN=$2; shift 2 ;;
    --gen-args) GEN_ARGS=$2; shift 2 ;;
    --compress) COMPRESS=$2; shift 2 ;;
    --deexcite) DEEX=geant4
                if [[ ${2:-} == geant4 || ${2:-} == beagle ]]; then DEEX=$2; shift; fi
                shift ;;
    --stage) STAGE=$2; shift 2 ;;
    -h|--help) sed -n '2,/^set -euo/p' "$0" | sed '$d'; exit 0 ;;
    *) echo "unknown option $1" >&2; exit 2 ;;
  esac
done
die() { echo "run_eic_sim: $*" >&2; exit 2; }
[[ -n $BEAMS && -n $ION && -n $OUT ]] || die "--beams, --ion and --out are required (see --help)"
[[ $STAGE == all || $STAGE == gen ]] && [[ -z $ATTEMPTS ]] && die "--attempts is required"
[[ $Q2 =~ ^[0-9.]+:[0-9.]+$ ]] && awk -F: '{exit !($1 >= 0 && $2 > $1)}' <<< "$Q2" || die "--q2 must be MIN:MAX with 0 <= MIN < MAX"
case "$ION" in C12) Z=6 N=6 ;; d) Z=1 N=1 ;; He4) Z=2 N=2 ;; *) die "--ion must be C12, d or He4 (Z/A = 1/2)" ;; esac
case "$BEAMS" in
  5x41)   : "${EPIC_CONFIG:=epic_craterlake_5x41}" "${FIELD=fields/beamline_5x41_He4.xml}" "${AB_PRESET:=ip6_hidiv_41x5}" ;;
  10x110) : "${EPIC_CONFIG:=epic_craterlake_10x110_H2}" "${FIELD=}" "${AB_PRESET:=ip6_eAu_110x10}" ;;
  *) die "--beams must be 5x41 or 10x110 (the light-ion lattices in epic)" ;;
esac
: "${EIC_SHELL:=$HOME/eic/eic-shell}" "${FIELDMAPS:=$HOME/eic/fieldmaps}"
[[ -n ${EPIC_INSTALL:-} ]] || die "set EPIC_INSTALL to an epic install prefix"
if [[ ! -f $EPIC_INSTALL/bin/thisepic.sh ]]; then
  # container-only installs (e.g. /opt/detector/epic-main) are checked inside eic-shell
  [[ $EPIC_INSTALL == /opt/* ]] || die "EPIC_INSTALL=$EPIC_INSTALL has no bin/thisepic.sh"
  echo "run_eic_sim: $EPIC_INSTALL is not visible on the host; assuming it exists in the container" >&2
fi
[[ -x $EIC_SHELL ]] || die "EIC_SHELL=$EIC_SHELL is not executable"
[[ -d $FIELDMAPS ]] || die "FIELDMAPS=$FIELDMAPS is not a directory"
if [[ $STAGE == all || $STAGE == reco || $STAGE == resume ]]; then
  [[ -n ${EICRECON_INSTALL:-} && -f $EICRECON_INSTALL/bin/eicrecon-this.sh ]] || die "set EICRECON_INSTALL"
fi
DEEX_BIN=$GCF/tools/deexcite/build/deexcite
: "${BEAGLESYS:=$HOME/src/BeAGLE}" "${FLUPRO:=$HOME/src/beagle-deps/fluka}" "${LHAPATH:=$HOME/src/beagle-deps/lhapdf5/share/lhapdf/PDFsets}"
if [[ $DEEX != 0 ]]; then
  [[ $GEN == genQE || $GEN == genDIS ]] || die "--deexcite needs the A-2 remnant in the event: genQE or genDIS"
  if [[ $DEEX == geant4 ]]; then
    [[ -x $DEEX_BIN ]] || die "--deexcite: build tools/deexcite inside eic-shell first (see --help)"
  else
    DEEX_BIN=$GCF/tools/beagle_deexcite.py
    [[ -x $BEAGLESYS/BeAGLE ]] || die "--deexcite beagle: no BeAGLE executable in BEAGLESYS=$BEAGLESYS"
    [[ -f $FLUPRO/data/nuclear.bin ]] || die "--deexcite beagle: no FLUKA data in FLUPRO=$FLUPRO"
    [[ -f $LHAPATH/cteq6ll.LHpdf ]] || die "--deexcite beagle: no LHAPDF 5 sets in LHAPATH=$LHAPATH"
  fi
fi
if [[ $STAGE == all ]]; then
  [[ -e $OUT ]] && die "$OUT exists; use a new directory or --stage to continue a run"
  mkdir -p "$OUT"
else
  [[ -d $OUT ]] || die "$OUT does not exist"
fi
OUT=$(cd "$OUT" && pwd)
cd "$OUT"
log() { echo "[$(date +%T)] $*" | tee -a run.log; }

# ---------------------------------------------------------------- generation (host)
if [[ $STAGE == all || $STAGE == gen ]]; then
  set +u  # thisroot.sh reads unset variables
  # shellcheck disable=SC1091
  source "$GCF/setup_env.sh" > /dev/null
  set -u
  E_IRF=$("$GCF/tools/gcf2hepmc.py" --beams "$BEAMS" --ion "$ION" --print-ebeam-irf)
  printf 'QSq %s %s\n' "${Q2%%:*}" "${Q2##*:}" > ps.txt
  case "$GEN" in genQE|genQE_FSI|genDIS|genQE_3N_FSI) ;; *) die "--generator must be genQE, genQE_FSI, genDIS or genQE_3N_FSI" ;; esac
  # E* of the A-2 remnant (-M) only exists for an SRC pair in a nucleus with A > 4: for the
  # deuteron it would give the empty remnant a mass, and for 4He it would turn the A-2 = 2 pn
  # remnant into a deuteron with mass m_d + E* (about 17 MeV), which does not exist; pp and nn are
  # unbound either way. E* = 0 is the edge of the 0-10 MeV band Pybus et al. 2020 (4He(e,e'pN) with
  # GCF, arXiv:2003.02318, doi:10.1016/j.physletb.2020.135429) used as a systematic for the 4He
  # A-2 system, which gcfNucleus mirrors with Estar_max = 0.010 under -r; that band is not applied
  # here. The mean field (-F) has no pair, and genQE_3N_FSI has no -M.
  ESTAR="-M"; { (( Z + N <= 4 )) || [[ $GEN == genQE_3N_FSI || " $GEN_ARGS " == *" -F "* ]]; } && ESTAR=""
  if [[ $GEN == genQE_FSI || $GEN == genQE_3N_FSI ]]; then
    F=$GCF/src/programs/genQE_FSI
    export GXMLPATH="$F/config/genie_override:$GENIE/config/G18_10a:$GENIE/config" GMSGCONF="$F/config/quiet_messenger.xml"
  fi
  CONV_ARGS=()
  case "$GEN" in
    genQE|genQE_FSI)
      CMD=("$GCF/src/programs/genQE_FSI/build/$GEN" "$Z" "$N" "$E_IRF" gen.root "$ATTEMPTS" -P ps.txt $ESTAR -S "$SEED") ;;
    genDIS)
      CMD=("$GCF/build/programs/genDIS/genDIS" "$Z" "$N" "$E_IRF" "$GCF/src/libraries/cross_sections/CTEQ/ct18nn.00.pds"
           gen.root "$ATTEMPTS" -P ps.txt $ESTAR -S "$SEED") ;;
    genQE_3N_FSI)
      # the rest-frame electron scatters by a fraction of a degree; the window scales as 1/E_IRF
      # (0.1-1.5 deg at 5x41 and 0.02-0.28 deg at 10x110 cover the Q2 = 3-25 GeV^2 phase space)
      TW=$(awk -v e="$E_IRF" 'BEGIN { printf "%.4f:%.4f", 0.1 * 440.33 / e, 1.5 * 440.33 / e }')
      CMD=("$GCF/src/programs/genQE_3N/build/genQE_3N_FSI" "$E_IRF" 1 gen.root "$ATTEMPTS" -A $((Z + N)) -Z "$Z"
           -q "$Q2" -t "$TW" -S "$SEED")
      CONV_ARGS=(--ebeam-irf "$E_IRF") ;;   # the 3N tree does not store q
  esac
  # shellcheck disable=SC2086
  log "${CMD[*]} $GEN_ARGS"
  # shellcheck disable=SC2086
  "${CMD[@]}" $GEN_ARGS > generator.log 2>&1 || die "$GEN failed (see generator.log)"
  printf '%s\n' "${CMD[*]:1} $GEN_ARGS" | sed "s#$GCF/##g" > generator.txt
  printf '%s\n' "$GEN $ESTAR $GEN_ARGS" > generator_short.txt
  # shellcheck disable=SC2086
  "$GCF/tools/gcf2hepmc.py" gen.root gen.hepmc3 --beams "$BEAMS" --ion "$ION" $CUTS "${CONV_ARGS[@]}" ${MAXEV:+-n $MAXEV} | tee -a run.log
  # the vendored GENIE tree is rebuilt in place, so judge cleanliness on everything else
  { git -C "$GCF" rev-parse HEAD
    echo "modified_outside_genie=$(git -C "$GCF" status --porcelain -- . ':!src/programs/genQE_FSI/Generator-R-3_06_02' | wc -l)"
  } > gcf_version.txt
fi

# ---------------------------------------------------------------- simulation and reconstruction (eic-shell)
case $STAGE in all|gen|sim|reco|resume|record) ;; *) die "--stage must be all, gen, sim, reco, resume or record" ;; esac
if [[ $STAGE != gen && $STAGE != record ]]; then
  cat > container_stage.sh <<EOF
#!/bin/bash
# generated by run_eic_sim.sh; runs inside eic-shell
set -uo pipefail
OUT="$OUT"
cd "\$OUT"
STAGE=$STAGE CHUNK=$CHUNK JOBS=$JOBS RJOBS=$RJOBS RTHREADS=$RTHREADS COMPRESS=$COMPRESS SEED=$SEED FIELDMAPS="$FIELDMAPS"
EPIC_INSTALL="$EPIC_INSTALL" EPIC_CONFIG="$EPIC_CONFIG" FIELD="$FIELD" AB_PRESET="$AB_PRESET"
EICRECON_INSTALL="${EICRECON_INSTALL:-}" DEEX=$DEEX DEEX_BIN="$DEEX_BIN"
export BEAGLESYS="$BEAGLESYS" FLUPRO="$FLUPRO" LHAPATH="$LHAPATH"
EOF
  cat >> container_stage.sh <<'EOF'
log() { echo "[$(date +%T)] $*" | tee -a run.log; }
[[ -f $EPIC_INSTALL/bin/thisepic.sh ]] || { log "no $EPIC_INSTALL/bin/thisepic.sh in the container"; exit 3; }
source "$EPIC_INSTALL/bin/thisepic.sh" "$EPIC_CONFIG"
[[ $(readlink -f "$DETECTOR_PATH") == $(readlink -f "$EPIC_INSTALL/share/epic") ]] || { log "DETECTOR_PATH mismatch: $DETECTOR_PATH"; exit 3; }

# run step $1 for every chunk, at most JOBS at a time, one work dir per chunk (the FileLoader
# writes calibrations/ into the cwd, and eicrecon reads the material map from there)
for_chunks() {
  local step=$1 max=$2
  for f in chunk_*.hepmc; do
    local i=${f//[!0-9]/}
    [[ $STAGE == resume && $(cat "exit_${step}_$i" 2>/dev/null) == 0 ]] && continue
    mkdir -p "work_$i"; ln -sfn "$FIELDMAPS" "work_$i/fieldmaps"
    while [ "$(jobs -rp | wc -l)" -ge "$max" ]; do wait -n || true; done
    { if (cd "work_$i" && "$step" "$i"); then echo 0; else echo 1; fi > "exit_${step}_$i"; } &
  done
  wait
  local n fail; n=$(ls exit_"${step}"_* 2>/dev/null | wc -l); fail=$(grep -LxF 0 exit_"${step}"_* 2>/dev/null | wc -l)
  log "$step: $n chunks, $fail failed"
  [ "$fail" -eq 0 ]
}
# Jobs write to local scratch and copy one recompressed file to the run directory (fewer, smaller
# writes to the shared disk). hadd -f<algo*100+level> rewrites every tree with the new compression.
LOCAL=$(mktemp -d "${TMPDIR:-/tmp}/run_eic_sim.XXXXXX")
trap 'rm -rf "$LOCAL"' EXIT
land() { # local_file final_path
  if [[ $COMPRESS != 0 ]]; then
    hadd -f"$COMPRESS" "$1.z" "$1" > "$1.hadd.log" 2>&1 || { cat "$1.hadd.log"; return 1; }
    mv -f "$1.z" "$2" && rm -f "$1" "$1.hadd.log"
  else
    mv -f "$1" "$2"
  fi
}
sim() {
  local out=$LOCAL/sim_$1.edm4hep.root
  npsim --compactFile ../geometry.xml --inputFiles "../chunk_$1.hepmc" --outputFile "$out" \
    --numberOfEvents -1 --random.seed $((SEED * 100000 + 10#$1 + 1)) > "../npsim_$1.log" 2>&1 &&
  land "$out" "../sim_$1.edm4hep.root"
}
RECO_FLAGS=(
  -PRPOTS:ForwardRomanPotStaticRecParticles:requireBeamProton=false
  -PRPOTS:ForwardRomanPotStaticRecParticles:requireMatchingMatrix=false
  -PRPOTS:ForwardRomanPotRecParticles:requireBeamProton=false
  -PFOFFMTRK:ForwardOffMRecParticles:requireBeamProton=false
  -PFOFFMTRK:ForwardOffMRecParticles:requireMatchingMatrix=false
  # match barrel ScFi and imaging clusters by eta/phi only: at the backward barrel end the imaging
  # layers see about a third of the shower, so the stock 50% energy check drops the ScFi cluster and
  # the electron has no cluster (5x41 Q2 3-5: 96 -> 201 of 300 identified; eic/EICrecon#2995, #190)
  -PBEMC:EcalBarrelClustersWithoutShapes:energyRelTolerance=-1
)
reco() {
  local out=$LOCAL/reco_$1.edm4eic.root
  eicrecon "../sim_$1.edm4hep.root" -Pnthreads="$RTHREADS" -Pdd4hep:xml_files=../geometry.xml "${RECO_FLAGS[@]}" \
    -Ppodio:output_file="$out" > "../eicrecon_$1.log" 2>&1 &&
  land "$out" "../reco_$1.edm4eic.root"
}

if [[ $STAGE == all || $STAGE == sim ]]; then
  cp "$DETECTOR_PATH/$EPIC_CONFIG.xml" geometry.xml
  if [[ -n $FIELD ]]; then
    n=$(grep -c 'fields/beamline_[^"]*\.xml' geometry.xml)
    [[ $n == 1 ]] || { log "expected one beamline field include in $EPIC_CONFIG.xml, found $n"; exit 3; }
    sed -i -E "s#fields/beamline_[^\"]+\.xml#$FIELD#" geometry.xml
  fi
  log "geometry: $EPIC_CONFIG with $(grep -o 'fields/beamline_[^"]*\.xml' geometry.xml)"
  AB_IN=gen.hepmc3
  if [[ $DEEX != 0 ]]; then
    # de-excitation of the A-2 remnant (Geant4 or BeAGLE); its products follow all other particles,
    # so the analysis's role order (e', lead, recoil) is unchanged
    DEEX_ARGS=(); [[ $DEEX == beagle ]] && DEEX_ARGS=(--jobs $(( JOBS < 16 ? JOBS : 16 )))
    "$DEEX_BIN" gen.hepmc3 gen_deex.hepmc3 --seed "$SEED" "${DEEX_ARGS[@]}" > deexcite.log 2>&1 || { log "deexcite ($DEEX) failed (deexcite.log)"; exit 4; }
    log "$(grep '^deexcite:' deexcite.log)"
    AB_IN=gen_deex.hepmc3
  fi
  abconv -p "$AB_PRESET" -f hepmc3 --plot-off -o ab "$AB_IN" > abconv.log 2>&1 || { log "abconv failed"; exit 4; }
  log "abconv $AB_PRESET: $(grep -c '^E ' ab.hepmc) events"
  rm -f chunk_*.hepmc
  awk -v per="$CHUNK" '
    /^HepMC::Asciiv3-END/ { next }
    !started && !/^E / { hdr = hdr $0 "\n"; next }
    /^E / { started = 1
            if (n % per == 0) {
              if (out) { print "HepMC::Asciiv3-END_EVENT_LISTING" > out; close(out) }
              out = sprintf("chunk_%04d.hepmc", n / per); printf "%s", hdr > out }
            n++ }
    { print > out }
    END { print "HepMC::Asciiv3-END_EVENT_LISTING" > out }' ab.hepmc
  # fetch the FileLoader inputs once, then let the jobs link them from this cache
  mkdir -p prefetch; ln -sfn "$FIELDMAPS" prefetch/fieldmaps
  (cd prefetch && npsim --compactFile ../geometry.xml -G --gun.particle geantino --gun.energy "1*GeV" -N 1 \
     --outputFile prefetch.edm4hep.root > ../prefetch.log 2>&1) || { log "prefetch failed"; exit 4; }
fi
export DETECTOR_CACHE="$OUT/prefetch/calibrations${DETECTOR_CACHE:+:$DETECTOR_CACHE}"
if [[ $STAGE == resume ]]; then
  [[ -f geometry.xml && -d prefetch/calibrations ]] && ls chunk_*.hepmc > /dev/null 2>&1 \
    || { log "resume: no split chunks, geometry or prefetch here; rerun with --stage sim"; exit 3; }
  log "resume: $(grep -lxF 0 exit_sim_* 2>/dev/null | wc -l) of $(ls chunk_*.hepmc | wc -l) chunks already simulated"
fi
if [[ $STAGE == all || $STAGE == sim || $STAGE == resume ]]; then
  log "npsim: $(ls chunk_*.hepmc | wc -l) chunks of $CHUNK, $JOBS parallel"
  for_chunks sim "$JOBS" || log "WARNING: npsim failures (see exit_sim_*)"
fi
if [[ $STAGE == all || $STAGE == reco || $STAGE == resume ]]; then
  [[ $STAGE == resume ]] && for f in chunk_*.hepmc; do  # re-reconstruct chunks whose sim was redone
    i=${f//[!0-9]/}; [[ sim_$i.edm4hep.root -nt exit_reco_$i ]] && rm -f "exit_reco_$i"; done
  source "$EICRECON_INSTALL/bin/eicrecon-this.sh"
  printf '%s\n' "${RECO_FLAGS[@]}" > reco_flags.txt
  log "eicrecon: $(command -v eicrecon), $RJOBS parallel, $RTHREADS threads/job"
  for_chunks reco "$RJOBS" || log "WARNING: eicrecon failures (see exit_reco_*)"
fi

# versions and inputs for run.json
mm=$(ls work_*/calibrations/materials-map.cbor 2>/dev/null | head -1)
{
  # the image (eic_xl nightly tags move) and the Spack versions it was built with
  echo "container_image=${APPTAINER_CONTAINER:-${SINGULARITY_CONTAINER:-}}"
  echo "container_from=$(sed -n 's/.*deffile\.from": "\([^"]*\)".*/\1/p' /.singularity.d/labels.json 2>/dev/null)"
  echo "container_build_date=$(sed -n 's/.*build-date": "\([^"]*\)".*/\1/p' /.singularity.d/labels.json 2>/dev/null)"
  echo "jug_info_sha256=$(sha256sum /etc/jug_info 2>/dev/null | cut -c1-64)"
  for p in dd4hep npsim geant4 acts podio edm4hep edm4eic; do
    echo "container_$p=$(grep -oE " $p@[^~+ ]+" /etc/jug_info 2>/dev/null | head -1 | cut -d@ -f2)"
  done
  echo "abconv=$(abconv -v 2>&1 | tail -1)"
  echo "geant4=$(geant4-config --version 2>/dev/null)"
  echo "root=$(root-config --version 2>/dev/null)"
  echo "detector_path=$DETECTOR_PATH"
  echo "geometry_sha256=$(sha256sum geometry.xml | cut -c1-64)"
  echo "material_map=${mm:+$(readlink -f "$mm")}"
  echo "material_map_sha256=${mm:+$(sha256sum "$(readlink -f "$mm")" | cut -c1-64)}"
  echo "eicrecon=$(command -v eicrecon 2>/dev/null)"
} > container_info.txt
EOF
  printf '%s\n' "bash $OUT/container_stage.sh" | "$EIC_SHELL" > container_stage.log 2>&1 || true
  tail -n 20 run.log
fi

# ---------------------------------------------------------------- run record
python3 - "$OUT" "$EPIC_INSTALL" "${EICRECON_INSTALL:-}" "$BEAMS" "$ION" "$SEED" "$CHUNK" "$EPIC_CONFIG" "$FIELD" "$AB_PRESET" "$EIC_SHELL" <<'PY'
import json, re, subprocess, sys
from datetime import datetime
from pathlib import Path
out, epic, eicrecon, beams, ion, seed, chunk, config, field, preset, shell = sys.argv[1:]
o = Path(out)
def git(path):
    src = Path(path).resolve().parent if path else None
    if not src or not (src / ".git").exists():
        return None
    r = lambda *a: subprocess.run(["git", "-C", str(src), *a], capture_output=True, text=True).stdout.strip()
    return {"dir": str(src), "commit": r("rev-parse", "HEAD"), "describe": r("describe", "--always", "--dirty"),
            "branch": r("rev-parse", "--abbrev-ref", "HEAD")}
sif = re.findall(r"/[^ ${}:]+\.sif", Path(shell).read_text())
sif = [s for s in sif if Path(s).exists()]
image = {"path": sif[0], "size": Path(sif[0]).stat().st_size,
         "mtime": datetime.fromtimestamp(Path(sif[0]).stat().st_mtime).isoformat(timespec="seconds")} if sif else None
info = {}
if (o / "container_info.txt").exists():
    info = dict(l.split("=", 1) for l in (o / "container_info.txt").read_text().splitlines() if "=" in l)
expected_mm = Path(epic).resolve().parent / "calibrations" / "materials-map.cbor"
exp_sha = None
if expected_mm.exists():
    exp_sha = subprocess.run(["sha256sum", str(expected_mm.resolve())], capture_output=True, text=True).stdout[:64]
count = lambda pat, val: sum((p.read_text().strip() == val) for p in o.glob(pat))
run = {
    "beams": beams, "ion": ion, "genQE_seed": int(seed), "chunk_size": int(chunk),
    "generator": (o / "generator_short.txt").read_text().split() if (o / "generator_short.txt").exists()
                 else (o / "generator.txt").read_text().split() if (o / "generator.txt").exists() else ["genQE", "-M"],
    "generator_command": (o / "generator.txt").read_text().strip() if (o / "generator_short.txt").exists() else None,
    "q2_range_GeV2": [float(x) for x in (o / "ps.txt").read_text().split()[1:3]] if (o / "ps.txt").exists() else None,
    "gcf": (o / "gcf_version.txt").read_text().split() if (o / "gcf_version.txt").exists() else None,
    "converter": json.loads((o / "gen.hepmc3.json").read_text()) if (o / "gen.hepmc3.json").exists() else None,
    "deexcitation": json.loads((o / "gen_deex.hepmc3.json").read_text()) if (o / "gen_deex.hepmc3.json").exists() else None,
    "epic": {"install": epic, "config": config, "field_swap": field or None, "git": git(epic)},
    "eicrecon": {"install": eicrecon or None, "git": git(eicrecon) if eicrecon else None,
                 "flags": (o / "reco_flags.txt").read_text().split() if (o / "reco_flags.txt").exists() else None},
    "abconv_preset": preset,
    "container_image": image,
    "container": info,
    "material_map_matches_epic_checkout": (info.get("material_map_sha256") == exp_sha) if exp_sha and info.get("material_map_sha256") else None,
    "jobs": {"sim_ok": count("exit_sim_*", "0"), "sim_failed": count("exit_sim_*", "1"),
             "reco_ok": count("exit_reco_*", "0"), "reco_failed": count("exit_reco_*", "1")},
}
(o / "run.json").write_text(json.dumps(run, indent=2) + "\n")
print(f"run.json: sim {run['jobs']['sim_ok']} ok / {run['jobs']['sim_failed']} failed, "
      f"reco {run['jobs']['reco_ok']} ok / {run['jobs']['reco_failed']} failed, "
      f"material map matches epic checkout: {run['material_map_matches_epic_checkout']}")
PY
