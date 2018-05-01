#pragma once
#include "Function.h"

class Constr1 :
  public Function 
{
  double getVal(const Vec &v) override
  {
    return v[0] * v[0] + v[1] - 1.0;
  }
  Vec getGrad(const Vec &v) override
  {
    Vec grad(DIM);
    grad[0] = 2.0 * v[0];
    grad[1] = 1.0;
    return grad;
  }
};