#include "SimplexSolver.h"
#include <iostream>
#include <string>
#include <iostream>
#include <iomanip>
#define EPS 1e-6
using namespace std;
#include "SepHyperplaneSolver.h"
#include "Constr1.h"
#include "Constr2.h"
#include "Constr3.h"
/*
 * DIM is defined as global defines in project settings
 */


int main(void)
{
  Vec vC(DIM);
  vC[0] = 0.0;
  vC[1] = 0.0;
  vC[2] = 1.0;
  double eps = 1e-3;

  Vec trueRes(DIM-1);
  trueRes[0] = 0.589754512301458;
  trueRes[1] = 0.6521896152200695;



  OmegaFuncs omF;
  Constr1 c1F;
  Constr2 c2F;
  Constr3 c3F;
  omF.push_back(&c1F);
  omF.push_back(&c2F);
  omF.push_back(&c3F);

  SepHyperplaneSolver sHySlvr(vC, eps);
  sHySlvr.CreateInitialConstraints(omF);
  Vec tempRes = sHySlvr.Solve();
  Vec res(DIM - 1);
  for (int i = 0; i < DIM - 1; ++i)
    res[i] = tempRes[i];
  Vec err = res - trueRes;
  cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
  cout << "x: " << fixed << setprecision(10) << res[0] << endl;
  cout << "y: " << fixed << setprecision(10) << res[1] << endl;
  cout << "Steps count: " << sHySlvr.GetIterCnt() << endl;
  cout << "Error: " << err.Norm() << endl;
 }