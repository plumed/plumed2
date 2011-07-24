#include <vector>
#include <map>

using namespace std;

class Grid  
{
 
protected:
 vector<double> grid_;
 vector<double> min_,max_,dx_;  
 vector<unsigned> nbin_;
 unsigned maxsize_;
 
public:
 Grid(vector<double> gmin, vector<double> gmax,
      vector<unsigned> nbin, bool doclear=true);
 Grid(vector<double> gmin,vector<double> gmax,
      vector<double> dx, bool doclear=true);

 vector<double> getMin() const;
 vector<double> getMax() const;
 vector<double> getSide() const;
 unsigned dimension() const;
 
 vector<unsigned> getIndices(unsigned index) const;
 vector<unsigned> getIndices(vector<double> x) const;
 unsigned getIndex(vector<unsigned> indices) const;
 unsigned getIndex(vector<double> x) const;
 vector<double> getPoint(unsigned index) const;
 vector<double> getPoint(vector<unsigned> indices) const;
 vector<double> getPoint(vector<double> x) const;
 
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

class GridWithSpline : public Grid
{
  vector< vector <double> > der_;
  vector<bool> pbc_;
  
  public:
   GridWithSpline(vector<double> gmin, vector<double> gmax, 
                  vector<unsigned> nbin, vector<bool> pbc);
   GridWithSpline(vector<double> gmin, vector<double> gmax, 
                  vector<double> dx, vector<bool> pbc);


  double getValue(unsigned index); 
  double getValue(vector<unsigned> indices);
  double getValue(vector<double> x);
  double getValue(vector<double> x, vector<double>& der);
  void setValue(unsigned index, double value);
  void setValue(unsigned index, double value, vector<double>& der);
  void setValue(vector<unsigned> indices, double value);
  void setValue(vector<unsigned> indices, double value, vector<double>& der);
  void addValue(unsigned index, double value);
  void addValue(unsigned index, double value, vector<double>& der);
  void addValue(vector<unsigned> indices, double value);
  void addValue(vector<unsigned> indices, double value, vector<double>& der);
  
  ~GridWithSpline(){};
};
  
class SparseGrid : public Grid
{

 protected:
 map<unsigned,double> map_;
 map<unsigned,double>::iterator it_;

 public:
 SparseGrid(vector<double> gmin, vector<double> gmax, vector<unsigned> nbin):
            Grid(gmin,gmax,nbin,false){};
 SparseGrid(vector<double> gmin,vector<double> gmax,vector<double> dx):
            Grid(gmin,gmax,dx,false){};

 unsigned size() const;
 void clear();
 double getRealSize() const;
 virtual double getValue(unsigned index); 
 virtual double getValue(vector<unsigned> indices);
 virtual double getValue(vector<double> x);
 virtual void setValue(unsigned index, double value);
 virtual void setValue(vector<unsigned> indices, double value);
 virtual void addValue(unsigned index, double value); 
 virtual void addValue(vector<unsigned> indices, double value);

 ~SparseGrid(){};
};

class SparseGridWithSpline : public SparseGrid
{
 map< unsigned,vector<double> > der_;
 map<unsigned,vector<double> >::iterator itder_;
 vector<bool> pbc_;

 public:
 SparseGridWithSpline(vector<double> gmin, vector<double> gmax,
            vector<unsigned> nbin, vector<bool> pbc):
            SparseGrid(gmin,gmax,nbin),pbc_(pbc){};
 SparseGridWithSpline(vector<double> gmin,vector<double> gmax,
            vector<double> dx, vector<bool> pbc):
            SparseGrid(gmin,gmax,dx),pbc_(pbc){};

 
  double getValue(unsigned index); 
  double getValue(vector<unsigned> indices);
  double getValue(vector<double> x);
  double getValue(vector<double> x, vector<double>& der);
  void setValue(unsigned index, double value);
  void setValue(unsigned index, double value, vector<double>& der);
  void setValue(vector<unsigned> indices, double value);
  void setValue(vector<unsigned> indices, double value, vector<double>& der);
  void addValue(unsigned index, double value);
  void addValue(unsigned index, double value, vector<double>& der);
  void addValue(vector<unsigned> indices, double value);
  void addValue(vector<unsigned> indices, double value, vector<double>& der);

 ~SparseGridWithSpline(){};
};

