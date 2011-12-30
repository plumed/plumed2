#ifndef __PLUMED_RMSD_h
#define __PLUMED_RMSD_h

#include "Vector.h"
#include "Log.h"
#include <vector>
#include <string>
#include "OptimalAlignment.h" 

namespace PLMD{

enum alignment_method_t {SIMPLE, OPTIMAL};

class PDB;

/// A class that implements RMSD calculations
class RMSD
{
  alignment_method_t alignment_method;
  std::vector<Vector> reference;
  std::vector<double> align;
  std::vector<double> displace;
  OptimalAlignment *myoptimalalignment;
  Log &log;
public:
/// initialize the log in the constructor
  RMSD(Log & log ): log(log){};
/// clear the structure
  void clear();
/// set reference, align and displace from input pdb structure
  void setFromPDB(const PDB&, std::string mytype);
/// set reference coordinates
  void setReference(const std::vector<Vector> & reference);
/// set weights
  void setAlign(const std::vector<double> & align);
/// set align
  void setDisplace(const std::vector<double> & displace);
/// 
  std::string getMethod();	
///
  double simpleAlignment(const  std::vector<double>  & align,
  		                     const  std::vector<double>  & displace,
  		                     const std::vector<Vector> & positions,
  		                     const std::vector<Vector> & reference ,
  		                     Log &log,
  		                     std::vector<Vector>  & derivatives);
/// Compute rmsd
  double calculate(const std::vector<Vector> & positions,std::vector<Vector> &derivatives);
};

}

#endif

