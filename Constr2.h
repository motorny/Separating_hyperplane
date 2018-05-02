#pragma once
#include "Function.h"

class Constr2 :
  public Function
{
  double getVal(const Vec &v) override
  {
    return v[0] + v[1]*v[1] / 3.0 - 1.0;
  }
  Vec getGrad(const Vec &v) override
  {
    Vec grad(DIM);
    grad[0] = 1.0;
    grad[1] = 2.0 * v[1] / 3.0;
    return grad;
  }
};
