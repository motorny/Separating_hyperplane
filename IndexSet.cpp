#include "IndexSet.h"
#include <iostream>
#include "Matrix.h"

IndexSet::IndexSet()
{
}

IndexSet::IndexSet(int cnt)
{
  for (int i = 0; i < cnt; ++i)
    data.push_back(i);
}

IndexSet::IndexSet(const IndexSet &is)
{
  data = is.data;
}

void IndexSet::Print(const char * s) const
{
  std::cout << s << std::endl;
  for (auto ind : data)
    std::cout << ind << " ";
  std::cout << "" << std::endl;
}

bool NextSet(int *a, int n, int m)
{
  for (int i = m - 1; i >= 0; --i)
    if (a[i] < n - m + i)
    {
      ++a[i];
      for (int j = i + 1; j < m; ++j)
        a[j] = a[j - 1] + 1;
      return true;
    }
  return false;
}
void IndexSet::CompleteToSize(const Matrix & m)
{
  int num = 0;
  IndexSet inv = (*this).GetInvertedSet(m.getColCnt());
  int invSetSize = inv.data.size();
  int *indArr = new int[invSetSize];

  int i = 0;
  for (auto ind : inv.data)
  {
    indArr[i++] = ind;
  }


  int toCompletecnt = m.getRowCnt() - data.size();

  int *indPosesInArr = new int[invSetSize];
  for (int i = 0; i < invSetSize; i++)
    indPosesInArr[i] = i;

  IndexSet part1 = IndexSet(*this);
  IndexSet part2 = IndexSet();
  IndexSet iSetFullRows = IndexSet(m.getRowCnt());
  do {
    part1.data.clear();
    part1 = IndexSet(*this);
    part2.data.clear();
    for (int i = 0; i < toCompletecnt; ++i)
    {
      part2.data.push_back(indArr[indPosesInArr[i]]);
    }


    part1.data.splice(part1.data.end(), part2.data);
    //part1.Print("Testing:");

    if (m.FullRang(iSetFullRows, part1))
    {
      for (int i = 0; i < toCompletecnt; ++i)
      {
        data.push_back(indArr[indPosesInArr[i]]);
      }
      delete indArr;
      delete indPosesInArr;
      return;
    }

  } while (NextSet(indPosesInArr, invSetSize, toCompletecnt));
}


void IndexSet::ChangeBasis(const Matrix & m, const IndexSet& iSetPositive, int * removedInd, int * addedInd, int * removedIndPos, Vec& d)
{
  IndexSet iSetOutter = (*this).GetInvertedSet(m.getColCnt());
  IndexSet iSetInner = IndexSet(*this);

  IndexSet iSetTest = IndexSet(*this);

  for (auto ind : iSetPositive.data)
  {
    iSetInner.data.remove(ind);
  }

  iSetInner.data.sort();
  iSetOutter.data.sort();

  iSetInner.Print("Inner");
  iSetOutter.Print("Outter");


  d.Print("WEIGHT:");
  iSetOutter.data.remove_if([d](int ind) { return d[ind] >= 0; });
   
  iSetOutter.Print("Outter");

  int iSetOutterSize = iSetOutter.data.size();
  int iSetInnerSize = iSetInner.data.size();
  int *indOutterArr = new int[iSetOutterSize];
  int *indInnerArr = new int[iSetInnerSize];

  int i = 0;
  for (auto ind : iSetOutter.data)
  {
    indOutterArr[i++] = ind;
  }

  int temp = indOutterArr[0];
  //for (int j = 0; j < iSetOutterSize; ++j)
  //{
  //  indOutterArr[j] = indOutterArr[(j + lastChangedPos + 1) % iSetOutterSize];
  //}
  //indOutterArr[iSetOutterSize - 1] = temp;

  i = 0;
  for (auto ind : iSetInner.data)
  {
    indInnerArr[i++] = ind;
  }

  int indPosO = 0;
  int indPosI = 0;
  int changedIndRelativePos = 0;

  IndexSet iSetFullRows = IndexSet(m.getRowCnt());

  while (indPosI < iSetOutterSize)
  {
    iSetTest.data.remove(indInnerArr[indPosI]);
    iSetTest.data.push_back(indOutterArr[indPosO]);
   // iSetTest.Print("Try:");
    if (m.FullRang(iSetFullRows, iSetTest))
    {
      std::list<int>::iterator iter = data.begin();
      for (; iter != data.end(); iter++)
      {
        if (*iter == indInnerArr[indPosI])
        {
          //*iter = indOutterArr[indPosO];
          break;
        }
        changedIndRelativePos++;
      }
      *removedInd = indInnerArr[indPosI];
      *addedInd = indOutterArr[indPosO];
      *removedIndPos = changedIndRelativePos;

      delete indOutterArr;
      delete indInnerArr;

      lastChangedPos++;
      return;
    }
    //revert changes
    iSetTest.data.remove(indOutterArr[indPosO]);
    iSetTest.data.push_back(indInnerArr[indPosI]);


    indPosO++;
    if (indPosO == iSetOutterSize)
    {
      indPosO = 0;
      indPosI++;
    }
  }
  std::cout << "Cant change basis" << std::endl;
}




IndexSet IndexSet::GetInvertedSet(int size_) const
{
  int num = 0;
  IndexSet temp(*this);
  temp.data.sort();
  std::list<int>::const_iterator it = temp.data.cbegin();
  IndexSet res;

  while (it != temp.data.end() && num < size_)
  {
    if (num < *it)
      res.data.push_back(num);
    else
      it++;
    num++;
  }

  if (it == temp.data.end())
    while (num < size_)
      res.data.push_back(num++);

  return res;
}
