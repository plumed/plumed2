#ifndef __PLUMED_Grid_h
#define __PLUMED_Grid_h

#include <vector>
#include <map>

namespace PLMD{ 

class Grid  
{
 std::vector<double> grid_;
 std::vector< std::vector<double> > der_;

protected:
 std::vector<double> min_,max_,dx_;  
 std::vector<unsigned> nbin_;
 std::vector<bool> pbc_;
 unsigned maxsize_, dimension_;
 bool dospline_, usederiv_;
 
 /// clear grid
 virtual void clear();
 
public:
 Grid(std::vector<double> gmin, std::vector<double> gmax, std::vector<unsigned> nbin, 
      std::vector<bool> pbc, bool dospline, bool usederiv, bool doclear=true);

/// get lower boundary
 std::vector<double> getMin() const;
/// get upper boundary
 std::vector<double> getMax() const;
/// get bin size
 std::vector<double> getDx() const;
/// get number of bins
 std::vector<unsigned> getNbin() const;
/// get grid dimension
 unsigned getDimension() const;
 
/// methods to handle grid indices 
 std::vector<unsigned> getIndices(unsigned index) const;
 std::vector<unsigned> getIndices(std::vector<double> x) const;
 unsigned getIndex(std::vector<unsigned> indices) const;
 unsigned getIndex(std::vector<double> x) const;
 std::vector<double> getPoint(unsigned index) const;
 std::vector<double> getPoint(std::vector<unsigned> indices) const;
 std::vector<double> getPoint(std::vector<double> x) const;

/// get neighbors
 std::vector<unsigned> getNeighbors(unsigned index,std::vector<unsigned> neigh);
 std::vector<unsigned> getNeighbors(std::vector<unsigned> indices,std::vector<unsigned> neigh);
 std::vector<unsigned> getNeighbors(std::vector<double> x,std::vector<unsigned> neigh);
/// get "neighbors" for spline
 std::vector<unsigned> getSplineNeighbors(std::vector<unsigned> indices);

/// get grid size
 virtual unsigned getSize() const;
/// get grid value
 virtual double getValue(unsigned index);
 virtual double getValue(std::vector<unsigned> indices);
 virtual double getValue(std::vector<double> x);
/// get grid value and derivatives
 virtual double getValueAndDerivatives(unsigned index, std::vector<double>& der); 
 virtual double getValueAndDerivatives(std::vector<unsigned> indices, std::vector<double>& der);
 virtual double getValueAndDerivatives(std::vector<double> x, std::vector<double>& der);

/// set grid value 
 virtual void setValue(unsigned index, double value);
 virtual void setValue(std::vector<unsigned> indices, double value);
/// set grid value and derivatives
 virtual void setValueAndDerivatives(unsigned index, double value, std::vector<double>& der);
 virtual void setValueAndDerivatives(std::vector<unsigned> indices, double value, std::vector<double>& der);
/// add to grid value
 virtual void addValue(unsigned index, double value); 
 virtual void addValue(std::vector<unsigned> indices, double value);
/// add to grid value and derivatives
 virtual void addValueAndDerivatives(unsigned index, double value, std::vector<double>& der); 
 virtual void addValueAndDerivatives(std::vector<unsigned> indices, double value, std::vector<double>& der); 

 virtual ~Grid(){};
};

  
class SparseGrid : public Grid
{

 std::map<unsigned,double> map_;
 std::map<unsigned,double>::iterator it_;
 std::map< unsigned,std::vector<double> > der_;
 std::map<unsigned,std::vector<double> >::iterator itder_;
 
 protected:
 void clear(); 
 
 public:
 SparseGrid(std::vector<double> gmin, std::vector<double> gmax, std::vector<unsigned> nbin,
            std::vector<bool> pbc, bool dospline, bool usederiv):
            Grid(gmin,gmax,nbin,pbc,dospline,usederiv,false){};
 
 unsigned getSize() const;
 double   getUsedSize() const;
 
 /// get grid value
 double getValue(unsigned index);
 double getValue(std::vector<unsigned> indices);
 double getValue(std::vector<double> x);
/// get grid value and derivatives
 double getValueAndDerivatives(unsigned index, std::vector<double>& der); 
 double getValueAndDerivatives(std::vector<unsigned> indices, std::vector<double>& der);
 double getValueAndDerivatives(std::vector<double> x, std::vector<double>& der);

/// set grid value 
 void setValue(unsigned index, double value);
 void setValue(std::vector<unsigned> indices, double value);
/// set grid value and derivatives
 void setValueAndDerivatives(unsigned index, double value, std::vector<double>& der);
 void setValueAndDerivatives(std::vector<unsigned> indices, double value, std::vector<double>& der);
/// add to grid value
 void addValue(unsigned index, double value); 
 void addValue(std::vector<unsigned> indices, double value);
/// add to grid value and derivatives
 void addValueAndDerivatives(unsigned index, double value, std::vector<double>& der); 
 void addValueAndDerivatives(std::vector<unsigned> indices, double value, std::vector<double>& der); 

 virtual ~SparseGrid(){};
};

}

#endif
