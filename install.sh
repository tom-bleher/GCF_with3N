#!/usr/bin/env bash
#
# install.sh — build the GCF_with3N generator suite, including GENIE FSI.
#
# Builds everything inside the checkout this script lives in. Nothing is
# installed into system locations and sudo is never used. Steps:
#   1. check prerequisites (ROOT, cmake, g++, make, ...; never installs them)
#   2. log4cpp 1.1.4 into a user prefix                  (skipped if present)
#   3. rebuild the vendored GENIE R-3_06_02 libraries in place for Linux
#      (the committed lib/ has no usable Linux libraries, only dangling
#      symlinks and dictionaries)                (skipped if valid libs exist)
#   4. generator suite without FSI         -> build/programs/<gen>/<gen>
#   5. 2N + FSI  (genQE_FSI project)       -> src/programs/genQE_FSI/build/
#   6. 3N + FSI  (genQE_3N project)        -> src/programs/genQE_3N/build/
#   7. write build/gcf_env.sh, the runtime environment for the FSI binaries
#   8. smoke test: generate a few events with each generator, check genT
#
# Usage:
#   ./install.sh [options]
#
# Options:
#   --prefix DIR   log4cpp install prefix                (default: $HOME/local)
#   -j, --jobs N   parallel build jobs                   (default: min(nproc, 32))
#   --root DIR     ROOT installation to use (sources DIR/bin/thisroot.sh)
#                  (default: root-config on PATH, else /opt/root-install)
#   --no-fsi       skip log4cpp and GENIE; build only the non-FSI generators
#   --force        clean-rebuild GENIE even if valid Linux libraries exist
#   --rebuild-log4cpp
#                  rebuild log4cpp even if it is installed (e.g. to replace
#                  an unoptimized build from before this script used -O2)
#   --no-test      skip the smoke test
#   -h, --help     show this help
#
# Per-step logs go to build/install-logs/. Re-running is cheap: finished
# steps are detected and skipped, and cmake only rebuilds what changed.
#
# NOTE: rebuilding GENIE rewrites ~130 *tracked* files under
# src/programs/genQE_FSI/Generator-R-3_06_02 (ROOT dictionaries, *.pcm,
# *.rootmap, Make.config, ...). They are machine-specific build products:
# do not commit them.
#
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$REPO_DIR/src"
BUILD_DIR="$REPO_DIR/build"              # suite build dir (gitignored)
FSI_DIR="$SRC_DIR/programs/genQE_FSI"
N3_DIR="$SRC_DIR/programs/genQE_3N"
GENIE_DIR="$FSI_DIR/Generator-R-3_06_02"
LOG_DIR="$BUILD_DIR/install-logs"
DEPS_DIR="$BUILD_DIR/deps"               # log4cpp tarball + source tree
SMOKE_DIR="$BUILD_DIR/smoke-test"
ENV_FILE="$BUILD_DIR/gcf_env.sh"

LOG4CPP_VERSION="1.1.4"
LOG4CPP_URL="https://sourceforge.net/projects/log4cpp/files/log4cpp-1.1.x%20%28new%29/log4cpp-1.1/log4cpp-${LOG4CPP_VERSION}.tar.gz/download"
LOG4CPP_SHA256="696113659e426540625274a8b251052cc04306d8ee5c42a0c7639f39ca90c9d6"

SUITE_BINS=(genQE/genQE genQE_3N/genQE_3N genDIS/genDIS genIon/genIon
            genPhoto/genPhoto genElectro/genElectro genCCQE/genCCQE
            sfSampler/sfSampler)
FSI2N_BINS=(genQE genQE_FSI SRC_analysis_2N)
FSI3N_BINS=(genQE_3N genQE_3N_FSI SRC_analysis_3N theta_plot theta_plot_fast
            theta_plot_sinelaw theta_scan)

# ----------------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------------
DIED=0
CUR_STEP="startup"
STEP=0
NSTEPS=0
T_START=$SECONDS

usage() {
    awk 'NR > 2 && /^#/ { sub(/^# ?/, ""); print; next } NR > 2 { exit }' "${BASH_SOURCE[0]}"
}
info() { printf '    %s\n' "$*"; }
warn() { printf 'WARNING: %s\n' "$*" >&2; }
die()  { DIED=1; printf '\nERROR: %s\n' "$*" >&2; exit 1; }

step() {
    STEP=$((STEP + 1))
    CUR_STEP="$*"
    STEP_T0=$SECONDS
    printf '\n==> [%d/%d] %s\n' "$STEP" "$NSTEPS" "$*"
}
step_done() { info "done ($((SECONDS - STEP_T0))s)"; }

on_exit() {
    local rc=$?
    if [[ $rc -ne 0 && $DIED -eq 0 ]]; then
        printf '\nERROR: install.sh failed (exit %d) during step "%s".\n' "$rc" "$CUR_STEP" >&2
        printf '       Step logs are in %s\n' "$LOG_DIR" >&2
    fi
}

# run_logged <log> <dir> <cmd...> — run cmd in dir, appending all output to
# log. On failure show the tail of the log and stop.
run_logged() {
    local log="$1" dir="$2" rc=0
    shift 2
    printf '\n+ (cd %s && %s)\n' "$dir" "$*" >>"$log"
    ( cd "$dir" && "$@" ) >>"$log" 2>&1 || rc=$?
    if [[ $rc -ne 0 ]]; then
        DIED=1
        printf '\nERROR: command failed (exit %d) in %s:\n  %s\n' "$rc" "$dir" "$*" >&2
        printf -- '---- last 25 lines of %s ----\n' "$log" >&2
        tail -n 25 "$log" >&2
        printf -- '----\nFull log: %s\n' "$log" >&2
        exit 1
    fi
}

# version_ge A B — true if version A >= B
version_ge() { [[ "$(printf '%s\n%s\n' "$2" "$1" | sort -V | head -n 1)" == "$2" ]]; }

# is_elf FILE — true if FILE (or the symlink's target) is an ELF binary
is_elf() { [[ -f "$1" && "$(head -c 4 "$1" 2>/dev/null)" == $'\177ELF' ]]; }

# source_root FILE — source thisroot.sh, dropping its harmless
# "ignored null byte" warning (thisroot.sh is not `set -u` clean)
source_root() {
    local err rc=0
    err="$(mktemp)"
    set +u
    # shellcheck disable=SC1090
    source "$1" 2>"$err" || rc=$?
    set -u
    grep -v 'ignored null byte' "$err" >&2 || true
    rm -f "$err"
    [[ $rc -eq 0 ]] || die "sourcing $1 failed (exit $rc)"
}

need_arg() { [[ $# -ge 2 && -n "$2" ]] || { echo "ERROR: $1 needs an argument (see --help)" >&2; exit 2; }; }

# ----------------------------------------------------------------------------
# Parse args
# ----------------------------------------------------------------------------
PREFIX="${HOME:?HOME is not set}/local"
JOBS=""
ROOT_OPT=""
DO_FSI=1
FORCE=0
REBUILD_LOG4CPP=0
DO_TEST=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --prefix)   need_arg "$@"; PREFIX="$2"; shift 2 ;;
        --prefix=*) PREFIX="${1#*=}"; shift ;;
        -j|--jobs)  need_arg "$@"; JOBS="$2"; shift 2 ;;
        -j*)        JOBS="${1#-j}"; shift ;;
        --jobs=*)   JOBS="${1#*=}"; shift ;;
        --root)     need_arg "$@"; ROOT_OPT="$2"; shift 2 ;;
        --root=*)   ROOT_OPT="${1#*=}"; shift ;;
        --no-fsi)   DO_FSI=0; shift ;;
        --force)    FORCE=1; shift ;;
        --rebuild-log4cpp) REBUILD_LOG4CPP=1; shift ;;
        --no-test)  DO_TEST=0; shift ;;
        -h|--help)  usage; exit 0 ;;
        *) echo "ERROR: unknown option '$1' (see --help)" >&2; exit 2 ;;
    esac
done

if [[ -z "$JOBS" ]]; then
    JOBS="$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)"
    if [[ "$JOBS" -gt 32 ]]; then JOBS=32; fi
fi
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --jobs needs a positive integer (got '$JOBS')" >&2; exit 2; }
trap on_exit EXIT

# Absolute, without trailing slash (it ends up in GENIE's Make.config).
[[ "$PREFIX" == /* ]] || PREFIX="$PWD/$PREFIX"
PREFIX="${PREFIX%/}"

NSTEPS=$(( 2 + (DO_FSI ? 5 : 0) + DO_TEST ))

echo "=================================================================="
echo " GCF_with3N install"
echo "   checkout      : $REPO_DIR"
echo "   FSI (GENIE)   : $([[ $DO_FSI -eq 1 ]] && echo "yes, log4cpp prefix $PREFIX" || echo "no (--no-fsi)")"
echo "   jobs          : $JOBS"
echo "   logs          : $LOG_DIR"
echo "=================================================================="

# ----------------------------------------------------------------------------
# 1. Prerequisites — report everything that is missing in one go.
# ----------------------------------------------------------------------------
step "Checking prerequisites"

[[ -f "$SRC_DIR/CMakeLists.txt" && -d "$GENIE_DIR/src" ]] \
    || die "$REPO_DIR does not look like a GCF_with3N checkout (no src/CMakeLists.txt or vendored GENIE)."
[[ "$(uname -s)" == "Linux" ]] \
    || die "this script supports Linux only (it rebuilds GENIE as ELF shared libraries)."
mkdir -p "$LOG_DIR"

problems=()
need() { command -v "$1" >/dev/null 2>&1 || problems+=("'$1' not found — $2"); }

# ROOT: --root, else whatever root-config is on PATH, else /opt/root-install.
THISROOT=""
if [[ -n "$ROOT_OPT" ]]; then
    THISROOT="${ROOT_OPT%/}/bin/thisroot.sh"
    [[ -f "$THISROOT" ]] || die "--root $ROOT_OPT: $THISROOT not found."
    source_root "$THISROOT"
elif command -v root-config >/dev/null 2>&1; then
    # already set up by the caller; only remember thisroot.sh for gcf_env.sh
    candidate="$(root-config --bindir)/thisroot.sh"
    if [[ -f "$candidate" ]]; then THISROOT="$candidate"; fi
elif [[ -f /opt/root-install/bin/thisroot.sh ]]; then
    THISROOT=/opt/root-install/bin/thisroot.sh
    source_root "$THISROOT"
fi
if command -v root-config >/dev/null 2>&1; then
    ROOT_VERSION="$(root-config --version)"
    [[ "$ROOT_VERSION" == 6.* ]] || problems+=("ROOT $ROOT_VERSION found, but ROOT 6 is required")
    if [[ $DO_FSI -eq 1 && "$(root-config --has-mathmore)" != "yes" ]]; then
        problems+=("ROOT was built without MathMore (GSL); GENIE needs it")
    fi
    info "ROOT     $ROOT_VERSION  ($(root-config --prefix))"
else
    problems+=("ROOT not found — no root-config on PATH and no /opt/root-install; pass --root DIR")
fi

if command -v cmake >/dev/null 2>&1; then
    out="$(cmake --version)"; out="${out%%$'\n'*}"; CMAKE_VERSION="${out##* }"
    version_ge "$CMAKE_VERSION" 3.6 || problems+=("cmake $CMAKE_VERSION is too old (need >= 3.6)")
    info "cmake    $CMAKE_VERSION"
else
    problems+=("'cmake' not found — need cmake >= 3.6")
fi
need g++ "need a C++17 compiler (GCC)"
command -v g++ >/dev/null 2>&1 && info "g++      $(g++ -dumpfullversion -dumpversion 2>/dev/null)"
if command -v make >/dev/null 2>&1; then
    out="$(make --version 2>/dev/null || true)"
    [[ "$out" == "GNU Make"* ]] || problems+=("'make' is not GNU make (GENIE's Makefiles need it)")
else
    problems+=("'make' not found — need GNU make")
fi

LOG4CPP_NEEDED=0
if [[ $DO_FSI -eq 1 ]]; then
    need perl "GENIE's configure script is written in perl"

    # libxml2 headers + library (GENIE dependency)
    XML2_INC="" XML2_LIBDIR=""
    inc_cands=() lib_cands=()
    if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists libxml-2.0 2>/dev/null; then
        for tok in $(pkg-config --cflags-only-I libxml-2.0); do inc_cands+=("${tok#-I}"); done
        lib_cands+=("$(pkg-config --variable=libdir libxml-2.0)")
    fi
    if command -v xml2-config >/dev/null 2>&1; then
        for tok in $(xml2-config --cflags); do [[ "$tok" == -I* ]] && inc_cands+=("${tok#-I}"); done
        for tok in $(xml2-config --libs);   do [[ "$tok" == -L* ]] && lib_cands+=("${tok#-L}"); done
    fi
    inc_cands+=(/usr/include/libxml2 /usr/local/include/libxml2)
    multiarch="$(gcc -print-multiarch 2>/dev/null || true)"
    [[ -n "$multiarch" ]] && lib_cands+=("/usr/lib/$multiarch")
    lib_cands+=(/usr/lib64 /usr/lib /usr/local/lib)
    for d in "${inc_cands[@]}"; do [[ -f "$d/libxml/tree.h" ]] && { XML2_INC="$d"; break; }; done
    for d in "${lib_cands[@]}"; do [[ -e "$d/libxml2.so" ]] && { XML2_LIBDIR="$d"; break; }; done
    if [[ -z "$XML2_INC" || -z "$XML2_LIBDIR" ]]; then
        problems+=("libxml2 development files not found (need libxml/tree.h and libxml2.so; e.g. Ubuntu package libxml2-dev)")
    else
        info "libxml2  $XML2_INC, $XML2_LIBDIR"
    fi

    if [[ $REBUILD_LOG4CPP -eq 0 && -f "$PREFIX/lib/liblog4cpp.so" && -f "$PREFIX/include/log4cpp/Category.hh" ]]; then
        info "log4cpp  already installed in $PREFIX (--rebuild-log4cpp to rebuild)"
    else
        LOG4CPP_NEEDED=1
        if [[ $REBUILD_LOG4CPP -eq 1 ]]; then
            info "log4cpp  will be rebuilt in $PREFIX (--rebuild-log4cpp)"
        else
            info "log4cpp  not in $PREFIX — will be built"
        fi
        command -v curl >/dev/null 2>&1 || command -v wget >/dev/null 2>&1 \
            || problems+=("neither 'curl' nor 'wget' found — needed to download log4cpp")
        need sha256sum "needed to verify the log4cpp download"
        need tar "needed to unpack log4cpp"
    fi
fi

if [[ ${#problems[@]} -gt 0 ]]; then
    printf '\nERROR: missing prerequisites:\n' >&2
    printf '  - %s\n' "${problems[@]}" >&2
    printf 'This script never installs system packages. Install the above (or ask your\nadmin), or point --root at a ROOT 6 installation, then re-run.\n' >&2
    DIED=1
    exit 1
fi
step_done

# Runtime env for the FSI binaries: GENIE's libraries find each other and
# log4cpp only through LD_LIBRARY_PATH (they carry no RPATH), so this is
# required at run time, not just at build time.
runtime_env_lines() {
    cat <<EOF
export GENIE="$GENIE_DIR"
export PATH="\$GENIE/bin:\$PATH"
export LD_LIBRARY_PATH="\$GENIE/lib:$PREFIX/lib\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}"
export GXMLPATH="$FSI_DIR/config/genie_override:\$GENIE/config/G18_10a:\$GENIE/config"
export GMSGCONF="$FSI_DIR/config/quiet_messenger.xml"
EOF
}

GENIE_VERSION="$(tr -d '[:space:]' <"$GENIE_DIR/VERSION")"
GENIE_NLIBS=0

# genie_libs_ok — every lib/libG*.so (the set the generators link against)
# must resolve to an ELF file; the committed tree only has dangling symlinks.
genie_libs_ok() {
    local lib="$GENIE_DIR/lib" f key n=0
    for key in libGFwMsg libGPhNuclSt libGPhHadTransp; do
        is_elf "$lib/$key-$GENIE_VERSION.so" || return 1
    done
    for f in "$lib"/libG*.so; do
        is_elf "$f" || return 1
        n=$((n + 1))
    done
    GENIE_NLIBS=$n
    [[ $n -gt 0 ]]
}

# genie_dirty_count — number of tracked files the GENIE build modified
genie_dirty_count() {
    if command -v git >/dev/null 2>&1 && git -C "$REPO_DIR" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
        git -C "$REPO_DIR" status --porcelain --untracked-files=no -- "$GENIE_DIR" | wc -l | tr -d ' '
    else
        echo 0
    fi
}

if [[ $DO_FSI -eq 1 ]]; then
    # ------------------------------------------------------------------------
    # 2. log4cpp into $PREFIX
    # ------------------------------------------------------------------------
    step "log4cpp $LOG4CPP_VERSION -> $PREFIX"
    if [[ $LOG4CPP_NEEDED -eq 0 ]]; then
        info "already installed, skipping"
    else
        log="$LOG_DIR/log4cpp.log"; : >"$log"
        info "log: $log"
        tarball="$DEPS_DIR/log4cpp-$LOG4CPP_VERSION.tar.gz"
        sha_ok() { [[ "$(sha256sum "$1" | cut -d ' ' -f 1)" == "$LOG4CPP_SHA256" ]]; }
        mkdir -p "$DEPS_DIR" "$PREFIX"
        if [[ -f "$tarball" ]] && sha_ok "$tarball"; then
            info "using cached $tarball"
        else
            info "downloading from SourceForge"
            rm -f "$tarball" "$tarball.part"
            if command -v curl >/dev/null 2>&1; then
                run_logged "$log" "$DEPS_DIR" curl -fsSL --retry 3 --connect-timeout 30 -o "$tarball.part" "$LOG4CPP_URL"
            else
                run_logged "$log" "$DEPS_DIR" wget --tries=3 -O "$tarball.part" "$LOG4CPP_URL"
            fi
            sha_ok "$tarball.part" || die "sha256 mismatch for the log4cpp download ($tarball.part):
  expected $LOG4CPP_SHA256
  got      $(sha256sum "$tarball.part" | cut -d ' ' -f 1)
The download is corrupt or SourceForge served an error page; delete it and retry."
            mv "$tarball.part" "$tarball"
            info "sha256 verified"
        fi
        # The tarball unpacks to ./log4cpp/; locate configure rather than guess.
        unpack="$DEPS_DIR/log4cpp-$LOG4CPP_VERSION"
        rm -rf "$unpack"; mkdir -p "$unpack"
        run_logged "$log" "$unpack" tar -xzf "$tarball"
        cfg="$(find "$unpack" -maxdepth 2 -type f -name configure -print -quit)"
        [[ -n "$cfg" ]] || die "no configure script found in $tarball"
        l4src="$(dirname "$cfg")"
        info "configure, make -j$JOBS, make install"
        # log4cpp's configure drops the usual default -O2, and GENIE calls into
        # log4cpp for every LOG statement, so an unoptimized build costs ~1.5x
        # in FSI generation (#22).
        run_logged "$log" "$l4src" ./configure --prefix="$PREFIX" --with-pthreads \
            CXXFLAGS=-O2 CFLAGS=-O2
        run_logged "$log" "$l4src" make -j"$JOBS"
        run_logged "$log" "$l4src" make install
        [[ -f "$PREFIX/lib/liblog4cpp.so" && -f "$PREFIX/include/log4cpp/Category.hh" ]] \
            || die "log4cpp installed, but $PREFIX/lib/liblog4cpp.so or include/log4cpp/Category.hh is missing; see $log"
    fi
    step_done

    # ------------------------------------------------------------------------
    # 3. GENIE R-3_06_02 libraries, rebuilt in place (no `make install`)
    # ------------------------------------------------------------------------
    step "GENIE $GENIE_VERSION libraries (in place, $GENIE_DIR)"
    export GENIE="$GENIE_DIR"
    export PATH="$GENIE/bin:$PATH"
    export LD_LIBRARY_PATH="$PREFIX/lib:$GENIE/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

    if [[ $FORCE -eq 0 ]] && genie_libs_ok; then
        info "$GENIE_NLIBS valid Linux libraries already in lib/, skipping (use --force to rebuild)"
    else
        log="$LOG_DIR/genie.log"; : >"$log"
        info "log: $log  (takes a few minutes)"
        if [[ $FORCE -eq 1 ]]; then
            info "--force: deleting GENIE object files and versioned libraries"
            find "$GENIE_DIR/src" -type f -name '*.o' -delete
            rm -f "$GENIE_DIR"/lib/libG*-"$GENIE_VERSION".so
        fi
        info "configure"
        run_logged "$log" "$GENIE_DIR" ./configure \
            --disable-pythia6 --disable-lhapdf5 --disable-flux-drivers --disable-geom-drivers \
            --with-compiler=gcc \
            --with-libxml2-inc="$XML2_INC" --with-libxml2-lib="$XML2_LIBDIR" \
            --with-log4cpp-inc="$PREFIX/include" --with-log4cpp-lib="$PREFIX/lib"

        # Plain `make` = `all` -> `apps`, whose prerequisites
        # (INITIAL_BUILD_TARGETS) build every library; `apps` itself then fails
        # in src/Apps (gNtpConv needs the flux drivers, gMakePhotonStrucFunc
        # needs LHAPDF), and we need no apps. So build exactly that target list,
        # read from GENIE's own Makefile, one top-level target at a time (they
        # are not ordered against each other); each sub-make still runs -j.
        targets="$(printf 'gcf-lib-targets:\n\t@echo "GCF_LIB_TARGETS: $(INITIAL_BUILD_TARGETS)"\n' \
                   | make -s --no-print-directory -C "$GENIE_DIR" -f Makefile -f - gcf-lib-targets 2>>"$log" \
                   | sed -n 's/^GCF_LIB_TARGETS: *//p')" || true
        [[ -n "$targets" ]] || die "could not read INITIAL_BUILD_TARGETS from $GENIE_DIR/Makefile; see $log"
        read -r -a target_list <<<"$targets"
        i=0
        for t in "${target_list[@]}"; do
            i=$((i + 1))
            info "$(printf '[%2d/%d] make -j%s %s' "$i" "${#target_list[@]}" "$JOBS" "$t")"
            run_logged "$log" "$GENIE_DIR" make -j"$JOBS" "$t"
        done
        genie_libs_ok || die "GENIE build finished but lib/ still has missing or non-ELF libG*.so; see $log"
        info "$GENIE_NLIBS libG*.so in lib/ (versioned libs + symlinks), all ELF"
        info "NOTE: this rebuild modified $(genie_dirty_count) tracked files under"
        info "      ${GENIE_DIR#"$REPO_DIR"/} — do NOT commit them."
    fi
    step_done
fi

# ----------------------------------------------------------------------------
# cmake helpers
# ----------------------------------------------------------------------------
# cmake_build <log> <source dir> <build dir> [cmake args...]
# (cd + `cmake <src>` rather than -S/-B, which needs cmake >= 3.13)
cmake_build() {
    local log="$1" src="$2" bld="$3"
    shift 3
    : >"$log"
    info "log: $log"
    mkdir -p "$bld"
    run_logged "$log" "$bld" cmake "$src" "$@"
    run_logged "$log" "$bld" cmake --build . -- -j"$JOBS"
}

# expect_bins <dir> <names...> — fail unless every binary was produced
expect_bins() {
    local dir="$1" b missing=()
    shift
    for b in "$@"; do [[ -x "$dir/$b" ]] || missing+=("$b"); done
    [[ ${#missing[@]} -eq 0 ]] || die "build finished but these binaries are missing in $dir: ${missing[*]}"
}

LOG4CPP_ARGS=(-DLOG4CPP_INCLUDE_DIR="$PREFIX/include" -DLOG4CPP_LIBRARY="$PREFIX/lib/liblog4cpp.so")

# ----------------------------------------------------------------------------
# 4. Generator suite (no FSI)
# ----------------------------------------------------------------------------
step "Generator suite, no FSI -> build/"
cmake_build "$LOG_DIR/cmake-suite.log" "$SRC_DIR" "$BUILD_DIR" -DUSE_FSI=OFF
expect_bins "$BUILD_DIR/programs" "${SUITE_BINS[@]}"
step_done

if [[ $DO_FSI -eq 1 ]]; then
    XML2_ARGS=(-DLIBXML2_INCLUDE_DIR="$XML2_INC" -DLIBXML2_LIBRARY="$XML2_LIBDIR/libxml2.so")

    # ------------------------------------------------------------------------
    # 5. 2N + FSI (USE_FSI defaults ON in this project and is fatal if unmet)
    # ------------------------------------------------------------------------
    step "genQE_FSI (2N + FSI) -> ${FSI_DIR#"$REPO_DIR"/}/build/"
    cmake_build "$LOG_DIR/cmake-genQE_FSI.log" "$FSI_DIR" "$FSI_DIR/build" \
        -DUSE_FSI=ON "${LOG4CPP_ARGS[@]}" "${XML2_ARGS[@]}"
    expect_bins "$FSI_DIR/build" "${FSI2N_BINS[@]}"
    step_done

    # ------------------------------------------------------------------------
    # 6. 3N + FSI. This CMakeLists does NOT fail when a dependency is missing:
    #    it prints a WARNING and quietly sets USE_FSI OFF, so catch that here.
    # ------------------------------------------------------------------------
    step "genQE_3N (3N + FSI) -> ${N3_DIR#"$REPO_DIR"/}/build/"
    log="$LOG_DIR/cmake-genQE_3N.log"
    cmake_build "$log" "$N3_DIR" "$N3_DIR/build" \
        -DUSE_FSI=ON "${LOG4CPP_ARGS[@]}" "${XML2_ARGS[@]}"
    # Two independent checks: the warning text (cmake word-wraps it, so squash
    # whitespace first) and whether the genQE_3N_FSI target was generated.
    if tr -s ' \n' '  ' <"$log" | grep -q 'disabling FSI build' \
       || ! grep -qs 'genQE_3N_FSI' "$N3_DIR/build/Makefile" "$N3_DIR/build/build.ninja"; then
        awk '/^CMake Warning/ { blk = $0; p = 1; next }
             p { blk = blk "\n" $0; if (/^$/) { if (blk ~ /disabling/) print blk; p = 0 } }' "$log" >&2
        die "genQE_3N's CMakeLists silently switched USE_FSI OFF, so genQE_3N_FSI was
not built. Check the GENIE libs/headers and log4cpp ($PREFIX); log: $log"
    fi
    expect_bins "$N3_DIR/build" "${FSI3N_BINS[@]}"
    step_done

    # ------------------------------------------------------------------------
    # 7. Runtime environment file (untracked; setup_env.sh is left untouched)
    # ------------------------------------------------------------------------
    step "Runtime environment -> ${ENV_FILE#"$REPO_DIR"/}"
    {
        echo "# GCF_with3N runtime environment, generated by install.sh on $(date '+%Y-%m-%d %H:%M')."
        echo "# Untracked (lives in build/); re-run ./install.sh to regenerate. Usage:"
        echo "#   source $ENV_FILE"
        if [[ -n "$THISROOT" ]]; then
            echo "source \"$THISROOT\""
        else
            echo "# ROOT: make sure $(command -v root-config) is on PATH"
        fi
        runtime_env_lines
    } >"$ENV_FILE"
    eval "$(runtime_env_lines)"
    info "wrote $ENV_FILE"
    if [[ "$PREFIX" != "$HOME/local" ]]; then
        info "(the tracked setup_env.sh assumes log4cpp in \$HOME/local, which is not your"
        info " --prefix: use build/gcf_env.sh instead)"
    fi
    step_done
fi

# ----------------------------------------------------------------------------
# 8. Smoke test. genQE_3N(_FSI) open array/he3_3Ndensity_*.dnst relative to
#    the CWD, so they must run from src/programs/genQE_3N.
# ----------------------------------------------------------------------------
# genT_entries <file.root> — print the entry count of tree genT (-1 if absent)
genT_entries() {
    local out
    out="$(root -l -b -q -e "TFile *f = TFile::Open(\"$1\"); TTree *t = f ? (TTree*)f->Get(\"genT\") : nullptr; printf(\"GCF_ENTRIES=%lld\\n\", t ? t->GetEntries() : -1LL);" 2>&1)" || true
    sed -n 's/^GCF_ENTRIES=//p' <<<"$out"
}

TIMEOUT=()
command -v timeout >/dev/null 2>&1 && TIMEOUT=(timeout 900)

# smoke <label> <run dir> <exact|atmost> <nEvents> <out.root> <cmd...>
#   exact : generator loops until nEvents weight>0 events (genQE_3N*)
#   atmost: nEvents attempts, only weight>0 ones stored    (genQE*)
smoke() {
    local label="$1" wd="$2" mode="$3" nev="$4" out="$5"
    shift 5
    local log="$LOG_DIR/smoke-$label.log" t0=$SECONDS attempt rc n
    : >"$log"
    rm -f "$out"
    for attempt in 1 2; do
        rc=0
        printf '+ (cd %s && %s)\n' "$wd" "$*" >>"$log"
        ( cd "$wd" && ${TIMEOUT[@]+"${TIMEOUT[@]}"} "$@" ) >>"$log" 2>&1 || rc=$?
        [[ $rc -eq 0 ]] && break
        # GENIE's cascade has a rare per-event segfault (run_3N_parallel.sh
        # retries for it too); a fresh auto-seed avoids it.
        if [[ $attempt -eq 1 && $rc -gt 128 && $rc -ne 124 ]]; then
            warn "$label died with signal $((rc - 128)); retrying once with a fresh seed"
        else
            die "smoke test $label failed (exit $rc$([[ $rc -eq 124 ]] && echo ', timed out')); see $log"
        fi
    done
    n="$(genT_entries "$out")"
    case "$mode" in
        exact)  [[ "$n" == "$nev" ]] ;;
        atmost) [[ "$n" =~ ^[0-9]+$ && "$n" -gt 0 && "$n" -le "$nev" ]] ;;
    esac || die "smoke test $label: genT has '${n:-?}' entries in $out (expected $([[ $mode == exact ]] && echo "$nev" || echo "1..$nev")); see $log"
    info "$(printf '%-13s %4s entries in genT (%s %s requested)  %3ds' "$label" "$n" \
        "$([[ $mode == exact ]] && echo "==" || echo "<=")" "$nev" $((SECONDS - t0)))"
}

if [[ $DO_TEST -eq 1 ]]; then
    step "Smoke test (outputs in ${SMOKE_DIR#"$REPO_DIR"/})"
    rm -rf "$SMOKE_DIR"; mkdir -p "$SMOKE_DIR"
    NEV=200
    smoke genQE "$SMOKE_DIR" atmost $NEV "$SMOKE_DIR/genQE.root" \
        "$BUILD_DIR/programs/genQE/genQE" 6 6 4.5 "$SMOKE_DIR/genQE.root" $NEV
    smoke genQE_3N "$N3_DIR" exact $NEV "$SMOKE_DIR/genQE_3N.root" \
        "$BUILD_DIR/programs/genQE_3N/genQE_3N" 6 1 "$SMOKE_DIR/genQE_3N.root" $NEV
    if [[ $DO_FSI -eq 1 ]]; then
        smoke genQE_FSI "$FSI_DIR" atmost $NEV "$SMOKE_DIR/genQE_FSI.root" \
            "$FSI_DIR/build/genQE_FSI" 6 6 4.5 "$SMOKE_DIR/genQE_FSI.root" $NEV
        smoke genQE_3N_FSI "$N3_DIR" exact $NEV "$SMOKE_DIR/genQE_3N_FSI.root" \
            "$N3_DIR/build/genQE_3N_FSI" 6 1 "$SMOKE_DIR/genQE_3N_FSI.root" $NEV
    fi
    step_done
fi

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
rel() { echo "${1#"$REPO_DIR"/}"; }
echo
echo "=================================================================="
echo " GCF_with3N build complete ($((SECONDS - T_START))s)   [$REPO_DIR]"
echo "------------------------------------------------------------------"
echo " Binaries (paths relative to the checkout):"
echo "   suite, no FSI : build/programs/<gen>/<gen>"
echo "                   ${SUITE_BINS[*]##*/}"
if [[ $DO_FSI -eq 1 ]]; then
    echo "   2N + FSI      : $(rel "$FSI_DIR")/build/"
    echo "                   ${FSI2N_BINS[*]}"
    echo "   3N + FSI      : $(rel "$N3_DIR")/build/"
    echo "                   ${FSI3N_BINS[*]}"
    echo
    echo " To run the FSI generators, in each new shell:"
    echo "   source $ENV_FILE"
    echo " which does exactly:"
    grep -v '^#' "$ENV_FILE" | sed 's/^/   /'
else
    echo
    echo " To run (no FSI), each new shell only needs ROOT:"
    if [[ -n "$THISROOT" ]]; then echo "   source $THISROOT"; else echo "   (root-config on PATH)"; fi
fi
echo
echo " Notes:"
echo "   * genQE_3N$([[ $DO_FSI -eq 1 ]] && echo " / genQE_3N_FSI") read array/he3_3Ndensity_*.dnst relative to"
echo "     the current directory: run them from $(rel "$N3_DIR"), e.g."
if [[ $DO_FSI -eq 1 ]]; then
    echo "       cd $N3_DIR && ./build/genQE_3N_FSI 6 1 /abs/path/out.root 1000"
else
    echo "       cd $N3_DIR && $BUILD_DIR/programs/genQE_3N/genQE_3N 6 1 /abs/path/out.root 1000"
fi
if [[ $DO_FSI -eq 1 ]]; then
    if [[ "$PREFIX" == "$HOME/local" ]]; then
        echo "   * The tracked setup_env.sh also works for this prefix, but does not set"
        echo "     GXMLPATH/GMSGCONF; build/gcf_env.sh is the complete version."
    else
        echo "   * The tracked setup_env.sh hard-codes log4cpp in \$HOME/local, which is"
        echo "     NOT your --prefix ($PREFIX): use build/gcf_env.sh."
    fi
    dirty="$(genie_dirty_count)"
    if [[ "$dirty" -gt 0 ]]; then
        echo "   * $dirty tracked files under $(rel "$GENIE_DIR") were modified by the"
        echo "     Linux GENIE build (dictionaries, *.pcm, *.rootmap, Make.config, ...)."
        echo "     They are machine-specific build products: do NOT commit them."
    fi
fi
echo "   * Logs: $LOG_DIR"
echo "=================================================================="
