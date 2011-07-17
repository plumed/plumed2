#include <vector>
#include <map>

using namespace std;

class Grid  
{
 vector<double> grid_;
 
protected:
 vector<double> min_,max_,dx_;  
 vector<unsigned> nbin_;
 unsigned maxsize_;
 bool dospline_;
 
public:
 Grid(vector<double> gmin,vector<double> gmax,
      vector<unsigned> nbin,bool dospline, bool doclear=true);

 double getMin(unsigned i) const;
 double getMax(unsigned i) const;
 double getSide(unsigned i) const;
 unsigned dimension() const;
 
 vector<unsigned> getIndices(unsigned index) const;
 vector<unsigned> getIndices(vector<double> x) const;
 unsigned getIndex(vector<unsigned> indices) const;
 unsigned getIndex(vector<double> x) const;
 vector<double> getPoint(unsigned index) const;
 vector<double> getPoint(vector<unsigned> indices) const;
 

//! a class that inherits from this, like SparseGrid, should override these methods
 virtual unsigned size() const;
 virtual void clear();
 virtual double getValue(unsigned index); 
 virtual double getValue(vector<unsigned> indices);
 virtual double getValue(vector<double> x);
 virtual void setValue(unsigned index, double value);
 virtual void setValue(vector<unsigned> indices, double value);
 virtual void addValue(unsigned index, double value); 
 virtual void addValue(vector<unsigned> indices, double value);

 ~Grid(){};
};

class SparseGrid : public Grid
{
 map<unsigned,double> map_;
 map<unsigned,double>::iterator it_;

 public:
 SparseGrid(vector<double> gmin,vector<double> gmax,vector<unsigned> nbin,bool dospline):
  Grid(gmin,gmax,nbin,dospline,false){};

 unsigned size() const;
 void clear();
 double getValue(unsigned index); 
 double getValue(vector<unsigned> indices);
 double getValue(vector<double> x);
 void setValue(unsigned index, double value);
 void setValue(vector<unsigned> indices, double value);
 void addValue(unsigned index, double value); 
 void addValue(vector<unsigned> indices, double value);

 ~SparseGrid(){};
};