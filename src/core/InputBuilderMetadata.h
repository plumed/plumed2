/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_core_InputBuilderMetadata_h
#define __PLUMED_core_InputBuilderMetadata_h

#include <iosfwd>
#include <map>
#include <string>
#include <vector>

namespace PLMD {

/// Write the registered action and command-line-tool metadata as JSON.
/// This is the shared serializer for gen_json and the in-memory metadata
/// commands. It enumerates the current registry on each call and emits UTF-8
/// JSON, including hidden keywords and explicit numbered-keyword metadata.
/// \param output Stream receiving the JSON document, including a final newline.
/// \param descriptions Optional action-summary overrides. If null, use the
/// embedded source-derived summaries. Missing entries give empty descriptions.
/// \param images Non-owning library handles used to resolve loaded registrations
/// in the same order as the action registry: later handles override earlier
/// handles and built-in registrations. An empty list selects built-in
/// registrations only. The caller must keep these libraries loaded throughout
/// the call; this function neither closes nor retains their handles.
/// \throws Exception If metadata is inconsistent or the output stream fails.
/// \throws std::ios_base::failure If the output stream has exceptions enabled.
void writeInputBuilderJSON(std::ostream& output,
                           const std::map<std::string,std::string>* descriptions=nullptr,
                           const std::vector<void*>& images= {});

} // namespace PLMD

#endif // __PLUMED_core_InputBuilderMetadata_h
