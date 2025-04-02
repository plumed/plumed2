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
#include "TrajectoryParser.h"

#include "FileBase.h"
#include "IFile.h"
#include "Tools.h"
#include "xdrfile/xdrfile_trr.h"
#include "xdrfile/xdrfile_xtc.h"


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

#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
static std::vector<molfile_plugin_t *> plugins;
static std::map <std::string, unsigned> pluginmap;
static int register_cb(void *v, vmdplugin_t *p) {
  //const char *key = p->name;
  const auto ret = pluginmap.insert ( std::pair<std::string,unsigned>(std::string(p->name),plugins.size()) );
  if (!ret.second) {
    //cerr<<"MOLFILE: found duplicate plugin for "<<key<<" : not inserted "<<endl;
  } else {
    //cerr<<"MOLFILE: loading plugin "<<key<<" number "<<plugins.size()-1<<endl;
    plugins.push_back(reinterpret_cast<molfile_plugin_t *>(p));
  }
  return VMDPLUGIN_SUCCESS;
}
#endif

void TrajectoryParser::registerKeywords(Keywords& keys) {
  keys.add("atoms","--ixyz","the trajectory in xyz format");
  keys.add("atoms","--igro","the trajectory in gro format");
  keys.add("atoms","--idlp4","the trajectory in DL_POLY_4 format");
  keys.add("atoms","--ixtc","the trajectory in xtc format (xdrfile implementation)");
  keys.add("atoms","--itrr","the trajectory in trr format (xdrfile implementation)");
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
  MOLFILE_INIT_ALL
  MOLFILE_REGISTER_ALL(NULL, register_cb)
  for(unsigned i=0; i<plugins.size(); i++) {
    std::string kk="--mf_"+std::string(plugins[i]->name);
    std::string mm=" molfile: the trajectory in "+std::string(plugins[i]->name)+" format " ;
    keys.add("atoms",kk,mm);
  }
#endif
}

std::optional<std::string> TrajectoryParser::init() {
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS

  molfile_timestep_t ts_in; // this is the structure that has the timestep
  // a std::vector<float> with the same scope as ts_in
  // it is necessary in order to store the pointer to ts_in.coords
  std::vector<float> ts_in_coords;
  ts_in.coords=ts_in_coords.data();
  ts_in.velocities=NULL;
  ts_in.A=-1; // we use this to check whether cell is provided or not
#endif

  int natoms;
  int command_line_natoms;
  std::string trajectoryFile;
  if(trajectory_fmt==trajfmt::molfile) {
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
    std::unique_ptr<std::lock_guard<std::mutex>> lck;
    if(api->is_reentrant==VMDPLUGIN_THREADUNSAFE) {
      lck=Tools::molfile_lock();
    }
    h_in = std::unique_ptr<void, molfile_deleter> {
      api->open_file_read(trajectoryFile.c_str(), trajectory_fmt_str.c_str(), &natoms),
      {api}
    };
    if(natoms==MOLFILE_NUMATOMS_UNKNOWN) {
      if(command_line_natoms>=0) {
        natoms=command_line_natoms;
      } else {
        return "this file format does not provide number of atoms; use --natoms on the command line";
      }
    }
    ts_in_coords.resize(3*natoms);
    ts_in.coords = ts_in_coords.data();
#endif
  } else if(trajectory_fmt==trajfmt::xdr_xtc || trajectory_fmt==trajfmt::xdr_trr) {
    xd.reset(xdrfile::xdrfile_open(trajectoryFile.c_str(),"r"));
    if(!xd) {
      return "ERROR: Error opening trajectory file "+trajectoryFile;
    }
    if(trajectory_fmt==trajfmt::xdr_xtc) {
      xdrfile::read_xtc_natoms(&trajectoryFile[0],&natoms);
    }
    if(trajectory_fmt==trajfmt::xdr_trr) {
      xdrfile::read_trr_natoms(&trajectoryFile[0],&natoms);
    }
  } else {
    fp.reset(std::fopen(trajectoryFile.c_str(),"r"));
    if(!fp) {
      return "ERROR: Error opening trajectory file "+trajectoryFile;
    }
  }


  return std::nullopt;
}

std::optional<std::string> TrajectoryParser::readFrame(
  int stride,
  bool dont_read_pbc,
  int & natoms,
  long long int &step,
  bool&localStep_toggle,
  real&localStep,
  real* masses,
  real* charges,
  real* coordinates,
  real* cell
) {
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS

  molfile_timestep_t ts_in; // this is the structure that has the timestep
  // a std::vector<float> with the same scope as ts_in
  // it is necessary in order to store the pointer to ts_in.coords
  std::vector<float> ts_in_coords;
  ts_in.coords=ts_in_coords.data();
  ts_in.velocities=NULL;
  ts_in.A=-1; // we use this to check whether cell is provided or not
#endif

  std::string line;

  //ensure that the file contains the next frame
  if(trajectory_fmt==trajfmt::molfile) {
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
    std::unique_ptr<std::lock_guard<std::mutex>> lck;
    if(api->is_reentrant==VMDPLUGIN_THREADUNSAFE) {
      lck=Tools::molfile_lock();
    }
    int rc;
    rc = api->read_next_timestep(h_in.get(), natoms, &ts_in);
    if(rc==MOLFILE_EOF) {
      return "Failure in getting the next frame with MOLFILE";
    }
#endif
  } else if(trajectory_fmt==trajfmt::xyz || trajectory_fmt==trajfmt::gro || trajectory_fmt==trajfmt::dlp4) {
    if(!Tools::getline(fp.get(),line)) {
      return "Failure in getting the next frame";
    }
  }
  //Here goes the `bool first_step=false;`
  if(trajectory_fmt!=trajfmt::molfile) {
    if(trajectory_fmt==trajfmt::xyz || trajectory_fmt==trajfmt::gro) {
      if(trajectory_fmt==trajfmt::gro)
        if(!Tools::getline(fp.get(),line)) {
          return "premature end of trajectory file";
        }
      std::sscanf(line.c_str(),"%100d",&natoms);
    }
    if(trajectory_fmt==trajfmt::dlp4) {
      char xa[9];
      int xb,xc,xd;
      double t;
      std::sscanf(line.c_str(),"%8s %lld %d %d %d %lf",xa,&step,&xb,&xc,&xd,&t);
      localStep_toggle=true;
      localStep=t;
      //   if (localstep) {
      //     p.cmd("setTimestep",real(t));
      //     localstep = false;
      //   }

    }
    ////////////////////////////////////////////////////////////////////////////
    // here the first step settings
    ////////////////////////////////////////////////////////////////////////////
    //now the actual reading of the trajectory
    if(trajectory_fmt==trajfmt::molfile) {
#ifdef __PLUMED_HAS_MOLFILE_PLUGINS
      if(!dont_read_pbc) {
        if(ts_in.A>0.0) { // this is negative if molfile does not provide box
          // info on the cell: convert using pbcset.tcl from pbctools in vmd distribution
          real cosBC=cos(real(ts_in.alpha)*pi/180.);
          //double sinBC=std::sin(ts_in.alpha*pi/180.);
          real cosAC=std::cos(real(ts_in.beta)*pi/180.);
          real cosAB=std::cos(real(ts_in.gamma)*pi/180.);
          real sinAB=std::sin(real(ts_in.gamma)*pi/180.);
          real Ax=real(ts_in.A);
          real Bx=real(ts_in.B)*cosAB;
          real By=real(ts_in.B)*sinAB;
          real Cx=real(ts_in.C)*cosAC;
          real Cy=(real(ts_in.C)*real(ts_in.B)*cosBC-Cx*Bx)/By;
          real Cz=std::sqrt(real(ts_in.C)*real(ts_in.C)-Cx*Cx-Cy*Cy);
          cell[0]=Ax/10.;
          cell[1]=0.;
          cell[2]=0.;
          cell[3]=Bx/10.;
          cell[4]=By/10.;
          cell[5]=0.;
          cell[6]=Cx/10.;
          cell[7]=Cy/10.;
          cell[8]=Cz/10.;
        } else {
          cell[0]=0.0;
          cell[1]=0.0;
          cell[2]=0.0;
          cell[3]=0.0;
          cell[4]=0.0;
          cell[5]=0.0;
          cell[6]=0.0;
          cell[7]=0.0;
          cell[8]=0.0;
        }
      }
      // info on coords
      // the order is xyzxyz...
      for(int i=0; i<3*natoms; i++) {
        coordinates[i]=real(ts_in.coords[i])/real(10.); //convert to nm
        //cerr<<"COOR "<<coordinates[i]<<endl;
      }
#endif
    } else if(trajectory_fmt==trajfmt::xdr_xtc || trajectory_fmt==trajfmt::xdr_trr) {
      int localstep;
      float time;
      xdrfile::matrix box;
      // here we cannot use a std::vector<rvec> since it does not compile.
      // we thus use a std::unique_ptr<rvec[]>
      auto pos=Tools::make_unique<xdrfile::rvec[]>(natoms);
      float prec,lambda;
      int ret=xdrfile::exdrOK;
      if(trajectory_fmt==trajfmt::xdr_xtc) {
        ret=xdrfile::read_xtc(xd.get(),natoms,&localstep,&time,box,pos.get(),&prec);
      }
      if(trajectory_fmt==trajfmt::xdr_trr) {
        ret=xdrfile::read_trr(xd.get(),natoms,&localstep,&time,&lambda,box,pos.get(),NULL,NULL);
      }
      if(stride==0) {
        step=localstep;
      }
      if(ret==xdrfile::exdrENDOFFILE) {
        return "EOF";
      }
      if(ret!=xdrfile::exdrOK) {
        return "EOF";
      }
      for(unsigned i=0; i<3; i++)
        for(unsigned j=0; j<3; j++) {
          cell[3*i+j]=box[i][j];
        }
      for(int i=0; i<natoms; i++)
        for(unsigned j=0; j<3; j++) {
          coordinates[3*i+j]=real(pos[i][j]);
        }
    } else {
      //headers
      if(trajectory_fmt==trajfmt::xyz) {
        if(!Tools::getline(fp.get(),line)) {
          return "premature end of trajectory file";
        }

        std::vector<double> celld(9,0.0);
        if(!dont_read_pbc) {
          std::vector<std::string> words;
          words=Tools::getWords(line);
          if(words.size()==3) {
            Tools::convert(words[0],celld[0]);
            Tools::convert(words[1],celld[4]);
            Tools::convert(words[2],celld[8]);
          } else if(words.size()==9) {
            Tools::convert(words[0],celld[0]);
            Tools::convert(words[1],celld[1]);
            Tools::convert(words[2],celld[2]);
            Tools::convert(words[3],celld[3]);
            Tools::convert(words[4],celld[4]);
            Tools::convert(words[5],celld[5]);
            Tools::convert(words[6],celld[6]);
            Tools::convert(words[7],celld[7]);
            Tools::convert(words[8],celld[8]);
          } else {
            return "needed box in second line of xyz file";
          }
        }
      }
      if(trajectory_fmt==trajfmt::dlp4) {
        std::vector<double> celld(9,0.0);
        if(!dont_read_pbc) {
          if(!Tools::getline(fp.get(),line)) {
            return "error reading vector a of cell";
          }
          std::sscanf(line.c_str(),"%lf %lf %lf",&celld[0],&celld[1],&celld[2]);
          if(!Tools::getline(fp.get(),line)) {
            return "error reading vector b of cell";
          }
          std::sscanf(line.c_str(),"%lf %lf %lf",&celld[3],&celld[4],&celld[5]);
          if(!Tools::getline(fp.get(),line)) {
            return "error reading vector c of cell";
          }
          std::sscanf(line.c_str(),"%lf %lf %lf",&celld[6],&celld[7],&celld[8]);
        }
      }
      int ddist=0;
      // Read coordinates
      for(int i=0; i<natoms; i++) {
        bool ok=Tools::getline(fp.get(),line);
        if(!ok) {
          return "premature end of trajectory file";
        }
        double cc[3];
        if(trajectory_fmt==trajfmt::xyz) {
          char dummy[1000];
          int ret=std::sscanf(line.c_str(),"%999s %100lf %100lf %100lf",dummy,&cc[0],&cc[1],&cc[2]);
          if(ret!=4) {
            return "cannot read line"+line;
          }
        } else if(trajectory_fmt==trajfmt::gro) {
          // do the gromacs way
          if(!i) {
            //
            // calculate the distance between dots (as in gromacs gmxlib/confio.c, routine get_w_conf )
            //
            const char      *p1, *p2, *p3;
            p1 = std::strchr(line.c_str(), '.');
            if (p1 == NULL) {
              return "seems there are no coordinates in the gro file";
            }
            p2 = std::strchr(&p1[1], '.');
            if (p2 == NULL) {
              return "seems there is only one coordinates in the gro file";
            }
            ddist = p2 - p1;
            p3 = std::strchr(&p2[1], '.');
            if (p3 == NULL) {
              return "seems there are only two coordinates in the gro file";
            }
            if (p3 - p2 != ddist) {
              return "not uniform spacing in fields in the gro file";
            }
          }
          Tools::convert(line.substr(20,ddist),cc[0]);
          Tools::convert(line.substr(20+ddist,ddist),cc[1]);
          Tools::convert(line.substr(20+ddist+ddist,ddist),cc[2]);
        } else if(trajectory_fmt==trajfmt::dlp4) {
          int lvl;
          char dummy[9];
          int idummy;
          double m,c;
          std::sscanf(line.c_str(),"%8s %d %lf %lf",dummy,&idummy,&m,&c);
          masses[i]=real(m);
          charges[i]=real(c);
          if(!Tools::getline(fp.get(),line)) {
            return "error reading coordinates";
          }
          std::sscanf(line.c_str(),"%lf %lf %lf",&cc[0],&cc[1],&cc[2]);
          cc[0]*=0.1;
          cc[1]*=0.1;
          cc[2]*=0.1;
          if(lvl>0) {
            if(!Tools::getline(fp.get(),line)) {
              return "error skipping velocities";
            }
          }
          if(lvl>1) {
            if(!Tools::getline(fp.get(),line)) {
              return "error skipping forces";
            }
          }
        } else {
          plumed_error();
        }
        // if(!debug_pd || ( i>=pd_start && i<pd_start+pd_nlocal) ) {
        //   coordinates[3*i]=real(cc[0]);
        //   coordinates[3*i+1]=real(cc[1]);
        //   coordinates[3*i+2]=real(cc[2]);
        // }
      }
      if(trajectory_fmt==trajfmt::gro) {
        if(!Tools::getline(fp.get(),line)) {
          return "premature end of trajectory file";
        }
        std::vector<std::string> words=Tools::getWords(line);
        if(words.size()<3) {
          return "cannot understand box format";
        }
        Tools::convert(words[0],cell[0]);
        Tools::convert(words[1],cell[4]);
        Tools::convert(words[2],cell[8]);
        if(words.size()>3) {
          Tools::convert(words[3],cell[1]);
        }
        if(words.size()>4) {
          Tools::convert(words[4],cell[2]);
        }
        if(words.size()>5) {
          Tools::convert(words[5],cell[3]);
        }
        if(words.size()>6) {
          Tools::convert(words[6],cell[5]);
        }
        if(words.size()>7) {
          Tools::convert(words[7],cell[6]);
        }
        if(words.size()>8) {
          Tools::convert(words[8],cell[7]);
        }
      }

    }
    return std::nullopt;
  }

}

} // namespace PLMD
