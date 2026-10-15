#ifndef __GCF_PROVENANCE_H__
#define __GCF_PROVENANCE_H__

// Provenance record written into every pipeline output file (as a TNamed
// named after the stage), so a number made from the file can be traced back
// to the exact command, seed and code version.

#include <ctime>
#include <sstream>
#include <string>

#ifndef GCF_GIT_VERSION
#define GCF_GIT_VERSION "unknown"
#endif

namespace gcfio {

inline std::string Provenance(const char *stage, int argc, char **argv,
                              unsigned long seed)
{
  std::ostringstream s;
  std::time_t now = std::time(nullptr);
  char when[64];
  std::strftime(when, sizeof when, "%Y-%m-%d %H:%M:%S %Z", std::localtime(&now));
  s << "stage: " << stage << "\n"
    << "version: " << GCF_GIT_VERSION << "\n"
    << "date: " << when << "\n"
    << "seed: " << seed << "\n"
    << "command:";
  for (int i = 0; i < argc; ++i) s << " " << argv[i];
  s << "\n";
  return s.str();
}

} // namespace gcfio

#endif
