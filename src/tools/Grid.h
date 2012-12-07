/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#ifndef __PLUMED_tools_Grid_h
#define __PLUMED_tools_Grid_h

#include <vector>
#include <string>
#include <map>

namespace PLMD{ 

class Value;
class PlumedIFile;
class PlumedOFile;
class KernelFunctions;

/// \ingroup TOOLBOX
class Grid  
{
 std::vector<double> grid_;
 std::vector< std::vector<double> > der_;
protected:
 std::string funcname;
 std::vector<std::string> argnames;
 std::vector<std::string> str_min_, str_max_;
 std::vector<double> min_,max_,dx_;  
 std::vector<unsigned> nbin_;
 std::vector<bool> pbc_;
 unsigned maxsize_, dimension_;
 bool dospline_, usederiv_;

 /// get "neighbors" for spline
 std::vector<unsigned> getSplineNeighbors(const std::vector<unsigned> & indices)const;

 /// clear grid
 virtual void clear();
 
public:
 Grid(const std::string& funcl, std::vector<Value*> args, const std::vector<std::string> & gmin, 
      const std::vector<std::string> & gmax, const std::vector<unsigned> & nbin, bool dospline, 
      bool usederiv, bool doclear=true);


/// get lower boundary
 std::vector<std::string> getMin() const;
/// get upper boundary
 std::vector<std::string> getMax() const;
/// get bin size
 std::vector<double> getDx() const;
/// get bin volume
 double getBinVolume() const;
/// get number of bins
 std::vector<unsigned> getNbin() const;
/// get if periodic
 std::vector<bool> getIsPeriodic() const;
/// get grid dimension
 unsigned getDimension() const;
 
/// methods to handle grid indices 
 std::vector<unsigned> getIndices(unsigned index) const;
 std::vector<unsigned> getIndices(const std::vector<double> & x) const;
 unsigned getIndex(const std::vector<unsigned> & indices) const;
 unsigned getIndex(const std::vector<double> & x) const;
 std::vector<double> getPoint(unsigned index) const;
 std::vector<double> getPoint(const std::vector<unsigned> & indices) const;
 std::vector<double> getPoint(const std::vector<double> & x) const;
/// faster versions relying on preallocated vectors
 void getPoint(unsigned index,std::vector<double> & point) const;
 void getPoint(const std::vector<unsigned> & indices,std::vector<double> & point) const;
 void getPoint(const std::vector<double> & x,std::vector<double> & point) const;

/// get neighbors
 std::vector<unsigned> getNeighbors(unsigned index,const std::vector<unsigned> & neigh) const;
 std::vector<unsigned> getNeighbors(const std::vector<unsigned> & indices,const std::vector<unsigned> & neigh) const;
 std::vector<unsigned> getNeighbors(const std::vector<double> & x,const std::vector<unsigned> & neigh) const;

/// write header for grid file
 void writeHeader(PlumedOFile& file);

/// read grid from file
 static Grid* create(const std::string&,std::vector<Value*>,PlumedIFile&,bool,bool,bool);
/// read grid from file and check boundaries are what is expected from input
 static Grid* create(const std::string&, std::vector<Value*>, PlumedIFile&,
                     const std::vector<std::string>&,const std::vector<std::string>&,
                     const std::vector<unsigned>&,bool,bool,bool); 
/// get grid size
 virtual unsigned getSize() const;
/// get grid value
 virtual double getValue(unsigned index) const;
 virtual double getValue(const std::vector<unsigned> & indices) const;
 virtual double getValue(const std::vector<double> & x) const;
/// get grid value and derivatives
 virtual double getValueAndDerivatives(unsigned index, std::vector<double>& der) const ;
 virtual double getValueAndDerivatives(const std::vector<unsigned> & indices, std::vector<double>& der) const;
 virtual double getValueAndDerivatives(const std::vector<double> & x, std::vector<double>& der) const;

/// set grid value 
 virtual void setValue(unsigned index, double value);
 virtual void setValue(const std::vector<unsigned> & indices, double value);
/// set grid value and derivatives
 virtual void setValueAndDerivatives(unsigned index, double value, std::vector<double>& der);
 virtual void setValueAndDerivatives(const std::vector<unsigned> & indices, double value, std::vector<double>& der);
/// add to grid value
 virtual void addValue(unsigned index, double value); 
 virtual void addValue(const std::vector<unsigned> & indices, double value);
/// add to grid value and derivatives
 virtual void addValueAndDerivatives(unsigned index, double value, std::vector<double>& der); 
 virtual void addValueAndDerivatives(const std::vector<unsigned> & indices, double value, std::vector<double>& der); 
/// Scale all grid values and derivatives by a constant factor
 virtual void scaleAllValuesAndDerivatives( const double& scalef );

/// add a kernel function to the grid
 void addKernel( const KernelFunctions& kernel );

/// dump grid on file
 virtual void writeToFile(PlumedOFile&);

 virtual ~Grid(){};
};

  
class SparseGrid : public Grid
{

 std::map<unsigned,double> map_;
 typedef std::map<unsigned,double>::const_iterator iterator;
 std::map< unsigned,std::vector<double> > der_;
 typedef std::map<unsigned,std::vector<double> >::const_iterator iterator_der;
 
 protected:
 void clear(); 
 
 public:
 SparseGrid(const std::string& funcl, std::vector<Value*> args, const std::vector<std::string> & gmin, 
            const std::vector<std::string> & gmax, 
            const std::vector<unsigned> & nbin, bool dospline, bool usederiv):
            Grid(funcl,args,gmin,gmax,nbin,dospline,usederiv,false){};
 
 unsigned getSize() const;
 unsigned getMaxSize() const;

/// this is to access to Grid:: version of these methods (allowing overloading of virtual methods)
 using Grid::getValue;
 using Grid::getValueAndDerivatives;
 using Grid::setValue;
 using Grid::setValueAndDerivatives;
 using Grid::addValue;
 using Grid::addValueAndDerivatives;
 
 /// get grid value
 double getValue(unsigned index) const;
/// get grid value and derivatives
 double getValueAndDerivatives(unsigned index, std::vector<double>& der) const; 

/// set grid value 
 void setValue(unsigned index, double value);
/// set grid value and derivatives
 void setValueAndDerivatives(unsigned index, double value, std::vector<double>& der);
/// add to grid value
 void addValue(unsigned index, double value); 
/// add to grid value and derivatives
 void addValueAndDerivatives(unsigned index, double value, std::vector<double>& der); 

/// dump grid on file
 void writeToFile(PlumedOFile&);

 virtual ~SparseGrid(){};
};

}

#endif
