#ifndef __PLUMED_Grid_h
#define __PLUMED_Grid_h

#include <vector>
#include <map>

using namespace std;

namespace PLMD{ 

class Grid  
{
 vector<double> grid_;
 vector< vector<double> > der_;

protected:
 vector<double> min_,max_,dx_;  
 vector<unsigned> nbin_;
 vector<bool> pbc_;
 unsigned maxsize_, dimension_;
 bool dospline_, usederiv_;

 /// get "neighbors" for spline
 vector<unsigned> getSplineNeighbors(vector<unsigned> indices);

 /// clear grid
 virtual void clear();
 
public:
 Grid(vector<double> gmin, vector<double> gmax, vector<unsigned> nbin, 
      vector<bool> pbc, bool dospline, bool usederiv, bool doclear=true);


/// get lower boundary
 vector<double> getMin() const;
/// get upper boundary
 vector<double> getMax() const;
/// get bin size
 vector<double> getDx() const;
/// get number of bins
 vector<unsigned> getNbin() const;
/// get grid dimension
 unsigned getDimension() const;
 
/// methods to handle grid indices 
 vector<unsigned> getIndices(unsigned index) const;
 vector<unsigned> getIndices(vector<double> x) const;
 unsigned getIndex(vector<unsigned> indices) const;
 unsigned getIndex(vector<double> x) const;
 vector<double> getPoint(unsigned index) const;
 vector<double> getPoint(vector<unsigned> indices) const;
 vector<double> getPoint(vector<double> x) const;

/// get neighbors
 vector<unsigned> getNeighbors(unsigned index,vector<unsigned> neigh);
 vector<unsigned> getNeighbors(vector<unsigned> indices,vector<unsigned> neigh);
 vector<unsigned> getNeighbors(vector<double> x,vector<unsigned> neigh);

/// get grid size
 virtual unsigned getSize() const;
/// get grid value
 virtual double getValue(unsigned index);
 virtual double getValue(vector<unsigned> indices);
 virtual double getValue(vector<double> x);
/// get grid value and derivatives
 virtual double getValueAndDerivatives(unsigned index, vector<double>& der); 
 virtual double getValueAndDerivatives(vector<unsigned> indices, vector<double>& der);
 virtual double getValueAndDerivatives(vector<double> x, vector<double>& der);

/// set grid value 
 virtual void setValue(unsigned index, double value);
 virtual void setValue(vector<unsigned> indices, double value);
/// set grid value and derivatives
 virtual void setValueAndDerivatives(unsigned index, double value, vector<double>& der);
 virtual void setValueAndDerivatives(vector<unsigned> indices, double value, vector<double>& der);
/// add to grid value
 virtual void addValue(unsigned index, double value); 
 virtual void addValue(vector<unsigned> indices, double value);
/// add to grid value and derivatives
 virtual void addValueAndDerivatives(unsigned index, double value, vector<double>& der); 
 virtual void addValueAndDerivatives(vector<unsigned> indices, double value, vector<double>& der); 

 virtual ~Grid(){};
};

  
class SparseGrid : public Grid
{

 map<unsigned,double> map_;
 map<unsigned,double>::iterator it_;
 map< unsigned,vector<double> > der_;
 map<unsigned,vector<double> >::iterator itder_;
 
 protected:
 void clear(); 
 
 public:
 SparseGrid(vector<double> gmin, vector<double> gmax, vector<unsigned> nbin,
            vector<bool> pbc, bool dospline, bool usederiv):
            Grid(gmin,gmax,nbin,pbc,dospline,usederiv,false){};
 
 unsigned getSize() const;
 double   getUsedSize() const;
 
 /// get grid value
 double getValue(unsigned index);
 double getValue(vector<unsigned> indices);
 double getValue(vector<double> x);
/// get grid value and derivatives
 double getValueAndDerivatives(unsigned index, vector<double>& der); 
 double getValueAndDerivatives(vector<unsigned> indices, vector<double>& der);
 double getValueAndDerivatives(vector<double> x, vector<double>& der);

/// set grid value 
 void setValue(unsigned index, double value);
 void setValue(vector<unsigned> indices, double value);
/// set grid value and derivatives
 void setValueAndDerivatives(unsigned index, double value, vector<double>& der);
 void setValueAndDerivatives(vector<unsigned> indices, double value, vector<double>& der);
/// add to grid value
 void addValue(unsigned index, double value); 
 void addValue(vector<unsigned> indices, double value);
/// add to grid value and derivatives
 void addValueAndDerivatives(unsigned index, double value, vector<double>& der); 
 void addValueAndDerivatives(vector<unsigned> indices, double value, vector<double>& der); 

 virtual ~SparseGrid(){};
};

}

#endif
