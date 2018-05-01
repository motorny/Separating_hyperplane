#pragma once
#include "Function.h"

class Constr3 :
  public Function
{
  double getVal(const Vec &v) override
  {
    return (v[0]-1)* (v[0] - 1) + (v[1] - 1)* (v[1] - 1) - v[2];
  }
  Vec getGrad(const Vec &v) override
  {
    Vec grad(DIM);
    grad[0] = 2.0 * (v[0] - 1);
    grad[1] = 2.0 * (v[1] - 1);
    grad[2] = -1.0;
    return grad;
  }
};
