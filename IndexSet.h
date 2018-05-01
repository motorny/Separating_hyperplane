#pragma once
#include <list>

class Matrix;
class Vec;

class IndexSet
{
public:
  std::list<int> data;
  int lastChangedPos = -1;
  IndexSet();
  /*
  Create indexSet filled with el-s from 0 to cnt-1
  e.g. 0, 1, 2, ....  , (cnt-1)
  */
  IndexSet(int cnt);
  IndexSet(const IndexSet &is);

  void Print(const char * s) const;

  void CompleteToSize(const Matrix &m);
  void ChangeBasis(const Matrix &m, const IndexSet& iSetPos, int * removedInd, int * addedInd, int * removedIndPos, Vec &d);
  IndexSet GetInvertedSet(int size_) const;

  int GetSize() const
  {
    return data.size();
  }
};

