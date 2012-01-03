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

 /// get "neighbors" for spline
 std::vector<unsigned> getSplineNeighbors(const std::vector<unsigned> & indices)const;

 /// clear grid
 virtual void clear();
 
public:
 Grid(const std::vector<double> & gmin, const std::vector<double> & gmax, const std::vector<unsigned> & nbin, 
      const std::vector<bool> & pbc, bool dospline, bool usederiv, bool doclear=true);


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
 std::vector<unsigned> getIndices(const std::vector<double> & x) const;
 unsigned getIndex(const std::vector<unsigned> & indices) const;
 unsigned getIndex(const std::vector<double> & x) const;
 std::vector<double> getPoint(unsigned index) const;
 std::vector<double> getPoint(const std::vector<unsigned> & indices) const;
 std::vector<double> getPoint(const std::vector<double> & x) const;

/// get neighbors
 std::vector<unsigned> getNeighbors(unsigned index,const std::vector<unsigned> & neigh) const;
 std::vector<unsigned> getNeighbors(const std::vector<unsigned> & indices,const std::vector<unsigned> & neigh) const;
 std::vector<unsigned> getNeighbors(const std::vector<double> & x,const std::vector<unsigned> & neigh) const;

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
 SparseGrid(const std::vector<double> & gmin, const std::vector<double> & gmax, const std::vector<unsigned> & nbin,
            const std::vector<bool> & pbc, bool dospline, bool usederiv):
            Grid(gmin,gmax,nbin,pbc,dospline,usederiv,false){};
 
 unsigned getSize() const;
 double   getUsedSize() const;
 
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

 virtual ~SparseGrid(){};
};

}

#endif
