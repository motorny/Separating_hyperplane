#pragma once

#include "Function.h"
#include <list>
class OmegaFuncs : public std::list<Function*> 
{
public:
  ;

  OmegaFuncs();
  ~OmegaFuncs();
  void Add(Function* f);
  double GetVal(const Vec& v) const;
  Vec GetGrad(const Vec& v) const;
};

