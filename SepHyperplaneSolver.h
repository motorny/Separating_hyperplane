#pragma once
#include "Matrix.h"
#include "SimplexSolver.h"
#include "Function.h"
#include "OmegaFuncs.h"



#define LOG_XK
#define LOG_YK
//#define LOG_AK
//#define LOG_BK
#define LOG_STEP_NORM



class SepHyperplaneSolver
{
private:
  Matrix mA;
  Vec vB, vC;
  OmegaFuncs omF;
  int iterCnt = 0;
  double eps;
public:
  SepHyperplaneSolver(const Vec& vC_, double eps_);
  void CreateInitialConstraints(const OmegaFuncs& omF_);
  Vec Solve();
  int GetIterCnt() const;
};

