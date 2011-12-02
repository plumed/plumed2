#ifndef __PLUMED_ColvarModifierFunctions_h
#define __PLUMED_ColvarModifierFunctions_h

#include <string>
#include <cassert>
#include <vector>
#include "SwitchingFunction.h"
#include "HistogramBead.h"
#include "ColvarModifier.h"

namespace PLMD {

class ColvarModifierMin : public ColvarModifier {
private:
   double beta;
   double total;
public:
   ColvarModifierMin(ColvarWithModifiers* act);
   void finishCalculation();
   double differential( const unsigned& ival, const double& value );
};

class ColvarModifierSum : public ColvarModifier {
private:
   double total;
public:
   ColvarModifierSum(ColvarWithModifiers* act);
   void finishCalculation();
   double differential( const unsigned& ival, const double& value );
};

class ColvarModifierMean : public ColvarModifier {
private:
   int ncv;
   double total;
public: 
   ColvarModifierMean (ColvarWithModifiers* act);
   void finishCalculation();
   double differential( const unsigned& ival, const double& value );
};

class ColvarModifierLess : public ColvarModifier {
private:
   SwitchingFunction sf;
   double total;
public:
   ColvarModifierLess(ColvarWithModifiers* act);
   void finishCalculation();
   double differential( const unsigned& ival, const double& value );
};

class ColvarModifierMore : public ColvarModifier {
private:
   SwitchingFunction sf;
   double total;
public: 
   ColvarModifierMore(ColvarWithModifiers* act);
   void finishCalculation();
   double differential( const unsigned& ival, const double& value );
};

class ColvarModifierHistogram : public ColvarModifier {
private:
   std::vector<HistogramBead> hf;
   std::vector<double> totals;
public:
   ColvarModifierHistogram(ColvarWithModifiers* act);
   void finishCalculation();
   double differential( const unsigned& ival, const double& value );
};



}

#endif
