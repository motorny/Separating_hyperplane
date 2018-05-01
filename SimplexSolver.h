#pragma once
#include "Vec.h"
#include "Matrix.h"

//#define NO_LOG
#define LOG_INITIAL_TASK
#define LOG_XK
#define LOG_BFS_INDS
#define LOG_B_MATRIX
#define LOG_DK
#define LOG_THETA
#define LOG_VU
#define LOG_IND_CHANGES

class SimplexSolver
{
public:
  static void ExtendTask(const Matrix& mA, const Vec& vB, Vec & vCE, Matrix & mAE, Vec & vBN, Vec & vX0E);
  static Vec Solve(const Vec& vC, const Matrix& mA, const Vec& vB, const Vec& vX0, Vec& DualSolution);
};

