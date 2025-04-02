/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2025 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#ifndef __PLUMED_tools_TrajectoryParser_h
#define __PLUMED_tools_TrajectoryParser_h
#include "Keywords.h"
#include "xdrfile/xdrfile.h"
#include "Tools.h"

#include <memory>
#include <optional>

// when using molfile plugin
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
#ifndef __PLUMED_HAS_EXTERNAL_MOLFILE_PLUGINS
/* Use the internal ones. Alternatively:
 *    ifeq (,$(findstring __PLUMED_HAS_EXTERNAL_MOLFILE_PLUGINS,$(CPPFLAGS)))
 *    CPPFLAGS+=-I../molfile
 */
#include "molfile/libmolfile_plugin.h"
#include "molfile/molfile_plugin.h"
using namespace PLMD::molfile;
#else
#include <libmolfile_plugin.h>
#include <molfile_plugin.h>
#endif
#endif


namespace PLMD {
// namespace molfile {
// struct molfile_plugin_t;
// } // namespace molfile



class TrajectoryParser {
  using real=double;
  struct file_deleter {
    void operator()(FILE* fp) {
      std::fclose(fp);
    }
  };

  struct xd_deleter {
    void operator()(PLMD::xdrfile::XDRFILE* xd) {
      xdrfile::xdrfile_close(xd);
    }
  };


  enum class trajfmt {
    molfile,
    xdr_xtc,
    xdr_trr,
    xyz,
    gro,
    dlp4
  } trajectory_fmt;
  std::string trajectory_fmt_str;



  std::unique_ptr<std::FILE, file_deleter> fp{nullptr,{}};
  std::unique_ptr<PLMD::xdrfile::XDRFILE, xd_deleter> xd{nullptr,{}};

#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
  molfile::molfile_plugin_t *api=NULL;
  struct molfile_deleter {
    molfile::molfile_plugin_t *api=NULL;
    void operator()(void* h_in) {
      if(h_in) {
        std::unique_ptr<std::lock_guard<std::mutex>> lck;
        if(api->is_reentrant==VMDPLUGIN_THREADUNSAFE) {
          lck=Tools::molfile_lock();
        }
        api->close_file_read(h_in);
      }
    }
  };
  std::unique_ptr<void, molfile_deleter> h_in{nullptr,{nullptr}};
#endif


  static void registerKeywords(Keywords& keys);
  std::optional<std::string> init();
  //if the opional is empty the frame is returned with success
  std::optional<std::string> readFrame(int stride,
                                       bool dont_read_pbc,
                                       int & natoms,
                                       long long int &step,
                                       bool&localStep_toggle,
                                       real&localStep,
                                       real* masses,
                                       real* charges,
                                       real* coordinates,
                                       real* cell );

//remember to add the
// else {			// from command line
//     celld=pbc_cli_box;
//   }
//   for(unsigned i=0; i<9; i++) {
//     cell[i]=real(celld[i]);
//   }

};
} //namespace PLMD
#endif //__PLUMED_tools_TrajectoryParser_h
