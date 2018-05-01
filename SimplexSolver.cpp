#include "SimplexSolver.h"
#include <iostream>

using namespace std;
void SimplexSolver::ExtendTask(const Matrix& mA, const Vec& vB, Vec & vCE, Matrix & mAE, Vec & vBN, Vec & vX0E)
{
  mAE = Matrix(mA.getRowCnt(), mA.getColCnt() + mA.getRowCnt());
  vCE = Vec(mA.getColCnt() + mA.getRowCnt());
  vBN = vB;
  vX0E = Vec(mA.getColCnt() + mA.getRowCnt());


  // copy main part of "A" matrix to new matrix
  for (int rI = 0; rI < mA.getRowCnt(); ++rI)
  {
    for (int cI = 0; cI < mA.getColCnt(); ++cI)
      mAE[rI][cI] = mA[rI][cI];
  }


  // invert rows of matrix, where corresponding component is negative
  for (int i = 0; i < mA.getRowCnt(); ++i)
    if (vBN[i] < 0)
    {
      vBN[i] = -vBN[i];
      for (int cI = 0; cI < mAE.getColCnt(); ++cI)
        mAE[i][cI] = -mAE[i][cI];
    }


  // append Identity matrix to "A" matrix
  for (int iDiag = 0; iDiag < mA.getRowCnt(); ++iDiag)
  {
    mAE[iDiag][mA.getColCnt() + iDiag] = 1.0;
  }

  // set initial basis
  for (int i = 0; i < mA.getRowCnt(); ++i)
    vX0E[mA.getColCnt() + i] = vBN[i];

  // set new goal vector
  for (int i = 0; i < mA.getRowCnt(); ++i)
    vCE[mA.getColCnt() + i] = 1.0;
}


Vec SimplexSolver::Solve(const Vec& vC, const Matrix& mA, const Vec& vB, const Vec& vX0, Vec& DualSolution)
{
  int stepNum = 0;

#ifdef LOG_INITIAL_TASK
  cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
  cout << "Initial task:" << endl;
  mA.Print("Matrix A:");
  vB.Print("Vector B:");
  vC.Print("Target vector C:");
  vX0.Print("Initial basic:");
  cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
#endif // LOG_INITIAL_TASK

  

  if (vC.getSize() != mA.getColCnt() || mA.getRowCnt() != vB.getSize())
  {
#ifndef NO_LOG
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "Incorrect Dimensions. Simplex Ended" << std::endl << std::endl;
#endif // !NO_LOG

    Vec err(1);
    err[0] = -1.0;
    return err;
  }
  IndexSet iSetBFS;
  IndexSet iSetFullRows(mA.getRowCnt());
  iSetBFS = vX0.GetPositiveInds();

#ifdef LOG_BFS_INDS
  iSetBFS.Print("Positive BFS (basic feasible solution) inds:");
#endif//LOG_BFS_INDS

  IndexSet iSetBFSPos = vX0.GetPositiveInds();
  bool degenterate = (int)iSetBFS.data.size() < mA.getRowCnt();

  iSetBFS.CompleteToSize(mA);

#ifdef LOG_BFS_INDS
  iSetBFS.Print("Completed BFS inds:");
#endif//LOG_BFS_INDS

  IndexSet iSetAntiBFS = iSetBFS.GetInvertedSet(mA.getColCnt());

#ifdef LOG_BFS_INDS
  iSetAntiBFS.Print("Lk (BFS inds complement):");
#endif//LOG_BFS_INDS

#ifdef LOG_B_MATRIX
  mA.Print("A matrix BFS:", iSetBFS);
#endif //LOG_B_MATRIX

  Matrix mB = mA.getInvertible(iSetFullRows, iSetBFS);

  Vec vXcur(vX0);
  Vec vXnext(vX0);
  Vec vd;

  Vec vU(mA.getColCnt());

  IndexSet iSetj;
  IndexSet iSetVDNeg;

  Matrix mTemp;

  for (;;)
  {
    stepNum++;
#if defined(LOG_XK) || defined(LOG_BFS_INDS) || defined(LOG_B_MATRIX) || defined(LOG_DK) || defined (LOG_THETA) || defined(LOG_VU) || defined(LOG_IND_CHANGES)
    std::cout << "___________________________________________" << std::endl;
    std::cout << "Starting step no: " << stepNum << std::endl;
#endif // 


    vXcur = vXnext;
#ifdef LOG_XK
    vXcur.Print("Curent BFS Vec: ");
#endif // LOG_XK


#ifdef LOG_B_MATRIX
    mB.Print("B matrix (Invertible):");
#endif //LOG_B_MATRIX

    vd = vC - vC.VecMulMatrix(iSetBFS, mB).VecMulMatrix(iSetFullRows, mA);

#ifdef LOG_DK
    vd.Print("dk vec:");
#endif // LOG_DK

    /*
     *Check for Optimum
     */
    iSetVDNeg = vd.ChoseNegative(iSetAntiBFS);
    if (iSetVDNeg.data.size() == 0)
    {
      //if all components are positeve then current vec is optimal

#ifndef NO_LOG
      std::cout << "-------------------------------------------" << std::endl;
      std::cout << "Solution found. Simplex Ended in " << stepNum - 1 << " steps" << std::endl << std::endl;
#endif // !NO_LOG

      DualSolution = vC.VecMulMatrix(iSetBFS, mB);
      return vXcur;
    }

    /*
     * Create and fill vU vector
     */
    iSetj.data.clear();
    int j = iSetVDNeg.data.front();
    iSetj.data.push_back(j);

    mTemp = mB.MatrixMulMatrix(iSetFullRows, iSetj, mA); // chose full rows set because mB is square m x m matirx

    //mTemp.Print("Temp Matix for u:");

    vU.SetZeros();
    int tempMInd = 0;
    for (auto ind : iSetBFS.data) {
      vU[ind] = mTemp[tempMInd++][0];
    }
    vU[j] = -1;

#ifdef LOG_VU
    vU.Print("vU vector:");
#endif // LOG_VU


    if (iSetBFS.data.size() == vU.ChoseNonPositive(iSetBFS).data.size())
    {
#ifndef NO_LOG
      std::cout << "-------------------------------------------" << std::endl;
      std::cout << "Unlimited. Simplex Ended" << std::endl << std::endl;
#endif // !NO_LOG

      Vec err(1);
      err[0] = -1.0;
      return err;
    }

    /*
     * Calculate Theta
     */
    int addedInd, removedInd;
    int changedIndRelativePos = 0;
    int minInd = -1;

    std::list<int>::iterator iterIndToChange;

    bool bfsChanged = false;

    IndexSet iSetTemp = vXcur.GetPositiveInds().GetInvertedSet(mA.getColCnt());//??????????????????????
    if (!degenterate || iSetTemp.data.size() == vU.ChoseNonPositive(iSetTemp).data.size())
    {
      IndexSet iSetPosvU = vU.ChosePositive(iSetBFS);
      minInd = iSetPosvU.data.front();
      double min = vXcur[minInd] / vU[minInd];
      for (auto ind : iSetPosvU.data)
      {
        if (vXcur[ind] / vU[ind] < min)
        {
          min = vXcur[ind] / vU[ind];
          minInd = ind;
        }
      }

#ifdef LOG_THETA
      std::cout << "Teta: " << min << std::endl;
#endif // LOG_THETA

      vXnext = vXnext - min * vU;

      std::list<int>::iterator iter = iSetBFS.data.begin();
      for (; iter != iSetBFS.data.end(); iter++)
      {
        if (*iter == minInd)
          break;
        changedIndRelativePos++;
      }
      iterIndToChange = iter;

      removedInd = minInd;
      addedInd = j;

    }
    else
    {
      std::cout << "$$$$$$$$ Changing Basis" << std::endl;
      iSetBFS.ChangeBasis(mA, iSetBFSPos, &removedInd, &addedInd, &changedIndRelativePos,vd);
      bfsChanged = true;
    }

#ifdef LOG_IND_CHANGES
    std::cout << "Removed Ind: " << removedInd << std::endl;
    std::cout << "Added Ind: " << addedInd << std::endl;
    std::cout << "Removed Ind pos: " << changedIndRelativePos << std::endl;
#endif // LOG_IND_CHANGES


    /*
     * Calculate new invertible
     */
    Matrix mF = Matrix(mA.getRowCnt(), mA.getRowCnt());
    for (int i = 0; i < mF.getRowCnt(); ++i)
      mF[i][i] = 1.0;

    int fI = 0;
    for (auto ind : iSetBFS.data)
      mF[fI++][changedIndRelativePos] = -vU[ind] / vU[removedInd];

    mF[changedIndRelativePos][changedIndRelativePos] = 1.0 / vU[removedInd];

    mB = mF * mB;

    if (!bfsChanged)
    {
      *iterIndToChange = j;
#ifdef LOG_BFS_INDS
      iSetBFS.Print("Next BFS:");
#endif//LOG_BFS_INDS
    }
    else
    {
      std::list<int>::iterator iter = iSetBFS.data.begin();
      for (; iter != iSetBFS.data.end(); iter++)
      {
        if (*iter == removedInd)
        {
          *iter = addedInd;
          break;
        }
      }
#ifdef LOG_BFS_INDS
      iSetBFS.Print("Next BFS:");
#endif//LOG_BFS_INDS
      mB = mA.getInvertible(iSetFullRows, iSetBFS);
    }
     //mA.Print("is invertible to:", iSetBFS);

    }// Main LOOP end

  }
