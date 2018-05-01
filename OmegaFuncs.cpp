#include "OmegaFuncs.h"
#include <iostream>
#include <algorithm>
using namespace std;

OmegaFuncs::OmegaFuncs()
{
}


OmegaFuncs::~OmegaFuncs()
{
}



double OmegaFuncs::GetVal(const Vec & v) const
{
  if (size() == 0)
  {
    cout << "Error! Functions list is empty!" << endl;
    return 0.0;
  }
  double max = front()->getVal(v);
  for (auto f : (*this))
  {
    if (f->getVal(v) > max)
      max = f->getVal(v);
  }
  return max;
}

Vec OmegaFuncs::GetGrad(const Vec & v) const
{
  if (size() == 0)
  {
    cout << "Error! Functions list is empty!" << endl;
    return Vec();
  }
  Function* fmaxP = front();
  double max = fmaxP->getVal(v);
  for (auto f : *this)
  {
    if (f->getVal(v) > max)
    {
      fmaxP = f;
      max = f->getVal(v);
    }
  }
  return fmaxP->getGrad(v);
}
