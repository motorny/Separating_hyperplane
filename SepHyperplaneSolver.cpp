#include "SepHyperplaneSolver.h"
#define EPS 1e-6
#include <iostream>
using namespace std;



SepHyperplaneSolver::SepHyperplaneSolver(const Vec& vC_, double eps_)
{
  eps = eps_;
  vC = vC_;
}

void SepHyperplaneSolver::CreateInitialConstraints(const OmegaFuncs & omF_)
{
  omF = omF_;
  mA = Matrix(3, 6);
  mA[0][0] = 1.0;
  mA[0][1] = -1.0;
  mA[1][2] = 1.0;
  mA[1][3] = -1.0;
  mA[2][4] = 1.0;
  mA[2][5] = -1.0;

  vB = Vec(6);
  vB[0] = 1.0;
  vB[1] = 2.0;
  vB[2] = 1.0;
  vB[3] = 3.0;
  vB[4] = 25.0;
  vB[5] = 0.0;
}

Vec SepHyperplaneSolver::Solve()
{
  Matrix mAE;
  Vec vBE, vCnewPositive, vY0E;
  Vec vXk;

  /*
  *Solve additional extended task
  */
  SimplexSolver::ExtendTask(mA, -1.0 * vC, vBE, mAE, vCnewPositive, vY0E);
  Vec BasicE = SimplexSolver::Solve(vBE, mAE, vCnewPositive, vY0E, vXk);

  if (BasicE[0] < 0)
  {
    Vec err(1);
    err[0] = -1.0;
    return err;
  }
  // Check gotten solution to analyze domain
  for (int i = mA.getColCnt(); i < BasicE.getSize(); ++i)
  {
    if (BasicE[i] > EPS)
    {
      cout << "Domain is empty!" << endl;
      Vec err(1);
      err[0] = -1.0;
      return err;
    }
  }

  // Form initial X vector for original task
  Vec Basic(mA.getColCnt());
  for (int i = 0; i < mA.getColCnt(); ++i)
    Basic[i] = BasicE[i];

  // Solve Original Task
  Vec vYk = SimplexSolver::Solve(vB, mA, -1.0 * vC, Basic, vXk);
  if (vYk[0] < 0)
  {
    Vec err(1);
    err[0] = -1.0;
    return err;
  }

#ifdef LOG_YK
  vYk.Print("Y on zero step:");
#endif // LOG_YK

#ifdef LOG_XK
  vXk.Print("X on zero step:");
#endif // LOG_XK

  Vec vLastXk;

  double stepEps = eps * 2.0;
  while (stepEps > eps)
  {
#if defined(LOG_XK) || defined(LOG_YK) || defined(LOG_B_MATRIX) || defined(LOG_DK) || defined (LOG_THETA) || defined(LOG_VU) || defined(LOG_IND_CHANGES)
    std::cout << "**************************************************************" << std::endl;
    std::cout << "Starting HYPERLANE step no: " << iterCnt<< std::endl;
#endif //
    iterCnt++;

    if (omF.GetVal(vXk) <= 0)
    {
      cout << "Point is in area. Finished in " << iterCnt - 1 << " steps!" << endl;
      return vXk;
    }

    // Starting New Step

    vLastXk = vXk;
    Vec vAk = omF.GetGrad(vXk);
#ifdef LOG_AK
    vAk.Print("New Constraint vector A:");
#endif // LOG_AK

    double bK = -omF.GetVal(vXk) + vAk * vXk;
#ifdef LOG_BK
    cout << "New Constraint b: " << bK << endl;
#endif // LOG_BK


    mA.Extend(0, 1);
    for (int rI = 0; rI < mA.getRowCnt(); ++rI)
    {
      mA[rI][mA.getColCnt() - 1] = vAk[rI];
    }

    vB.Extend(1);
    vB[vB.getSize() - 1] = bK;

    vYk.Extend(1);

    vYk = SimplexSolver::Solve(vB, mA, -1.0 * vC, vYk, vXk);
    
#ifdef LOG_YK
    vYk.Print("Yk Vector:");
#endif // LOG_YK

#ifdef LOG_XK
    vXk.Print("Xk Vector:");
#endif // LOG_XK

    Vec step = vXk - vLastXk;
    stepEps = step.Norm();
    cout << "Step norm: " << stepEps << endl;
  }
 
  return vXk;
}

int SepHyperplaneSolver::GetIterCnt() const
{
  return iterCnt;
}
