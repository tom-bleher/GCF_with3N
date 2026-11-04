# Writes gcf_version.h with the version of the source tree being built. Run at
# every build, not when CMake configures, so that an incremental build can't
# record a stale version in its output files. The file is only rewritten when
# the version changes, so an unchanged tree recompiles nothing.
#
#   cmake -DSRC_ROOT=<repo>/src -DOUT=<build>/gcf_version.h -P WriteVersion.cmake

execute_process(COMMAND git -C "${SRC_ROOT}" describe --always --tags
                OUTPUT_VARIABLE version OUTPUT_STRIP_TRAILING_WHITESPACE
                RESULT_VARIABLE describe_failed ERROR_QUIET)
if(describe_failed OR NOT version)
    set(version "unknown")
else()
    # Building the vendored GENIE modifies files that are tracked in git, so
    # they would mark every tree as dirty; only the code around them counts.
    execute_process(COMMAND git -C "${SRC_ROOT}" diff --quiet HEAD -- .
                            ":(exclude)programs/genQE_FSI/Generator-R-3_06_02"
                    RESULT_VARIABLE dirty ERROR_QUIET)
    if(dirty)
        set(version "${version}-dirty")
    endif()
endif()

set(content "#define GCF_GIT_VERSION \"${version}\"\n")
set(old "")
if(EXISTS "${OUT}")
    file(READ "${OUT}" old)
endif()
if(NOT content STREQUAL old)
    file(WRITE "${OUT}" "${content}")
endif()
