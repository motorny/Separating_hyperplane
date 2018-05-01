#include "Vec.h"
#include "Matrix.h"
using namespace std;
#define EPS 1e-6
#include <iostream>
#include <iomanip>
Vec::Vec()
{
}
Vec::Vec(int size_)
{
  data = new mT[size_];
  for (int i = 0; i < size_; ++i) {
    data[i] = 0.0;
  }
  size = size_;
}

Vec::Vec(const Vec & v)
{
  size = v.getSize();
  data = new mT[size];
  for (int i = 0; i < size; ++i)
    data[i] = v[i];
}

Vec & Vec::operator=(const Vec & v)
{
  if (size != v.getSize())
  {
    delete[] data;
    size = v.getSize();
    data = new mT[size];
  }

  for (int i = 0; i < size; ++i)
    data[i] = v[i];
  return (*this);
}

Vec::~Vec()
{
  delete[] data;
}




IndexSet Vec::GetPositiveInds() const
{
  IndexSet res;
  for (int i = 0; i < size; ++i)
  {
    if (data[i] > EPS)
      res.data.push_back(i);
  }
  return res;
}

Vec  Vec::VecMulMatrix(const IndexSet & iSetV, const Matrix & m)const
{
  if (iSetV.GetSize() != m.getRowCnt())
  {
    std::cout << "VEC * MATRIX DIMENSION ERROR" << std::endl;
    return Vec();
  }
  
  int resColCnt = m.getColCnt();
  Vec res(resColCnt);

  for (int i = 0; i < resColCnt; ++i)
  {
    mT sum = 0;
    int mInd = 0;
    for (auto ind : iSetV.data)
    {
      sum += data[ind] * m[mInd++][i];
    }
    res[i] = sum;
  }
  return res;
}

void Vec::SetZeros()
{
  for (int i = 0; i < size; ++i)
    data[i] = 0.0;
}

void Vec::Extend(int addRowCnt)
{
  int newSize = size + addRowCnt;
  mT* tempData = new mT[newSize];
  for (int i = 0; i < size; ++i)
    tempData[i] = data[i];
  for (int i = size; i < newSize; ++i)
    tempData[i] = 0.0;

  delete[] data;
  data = tempData;
  size = newSize;
}



void Vec::Print(const char * s) const
{
  std::cout << s << std::endl;
  for (int iRow = 0; iRow < size; iRow++)
  {
    std::cout << std::setw(5) << data[iRow] << endl;
  }
}

void Vec::Print(const char * s, const IndexSet & rowSet) const
{
  std::cout << s << std::endl;
  for (auto iRow : rowSet.data)
  {
    std::cout << std::setw(5) << data[iRow] << endl;
  }
}

IndexSet  Vec::ChoseNegative(const IndexSet & iSetV) const
{
  IndexSet res;
  for (auto ind : iSetV.data)
  {
    if (data[ind] < -EPS )
      res.data.push_back(ind);
  }
  return res;
}

IndexSet  Vec::ChosePositive(const IndexSet & iSetV) const
{
  IndexSet res;
  for (auto ind : iSetV.data)
  {
    if (data[ind] > 0)
      res.data.push_back(ind);
  }
  return res;
}

IndexSet  Vec::ChoseNonPositive(const IndexSet & iSetV) const
{
  IndexSet res;
  for (auto ind : iSetV.data)
  {
    if (data[ind] <= 0)
      res.data.push_back(ind);
  }
  return res;
}

mT Vec::Norm() const
{
   return sqrt((*this)*(*this));
}

mT & Vec::operator[](size_t ind)
{
  return data[ind];
}

mT Vec::operator[](size_t ind) const
{
  return data[ind];
}

Vec operator+(const Vec & v1, const Vec & v2)
{
  if(v1.getSize()!= v2.getSize())
  { 
    std::cout << "VEC DIMENSION ERROR" << std::endl;
    return Vec();
  }

  Vec res(v1.getSize());
  for (int i = 0; i < v1.getSize(); ++i)
    res[i] = v1[i] + v2[i];

  return res;
}

Vec operator-(const Vec & v1, const Vec & v2)
{
  if (v1.getSize() != v2.getSize())
  {
    std::cout << "VEC DIMENSION ERROR" << std::endl;
    return Vec();
  }

  Vec res(v1.getSize());
  for (int i = 0; i < v1.getSize(); ++i)
    res[i] = v1[i] - v2[i];

  return res;
}

Vec operator*(const Vec & v1, mT a)
{ 
  Vec res(v1.getSize());
  for (int i = 0; i < v1.getSize(); ++i)
    res[i] = v1[i] * a;
  return res;
}

Vec operator*(mT a, const Vec & v1)
{
  Vec res(v1.getSize());
  for (int i = 0; i < v1.getSize(); ++i)
    res[i] = v1[i] * a;
  return res;
}

mT operator*(const Vec & v1, const Vec & v2)
{
  mT sum = 0.0;
  if (v1.getSize() != v2.getSize())
  {
    std::cout << "VEC DIMENSION ERROR" << std::endl;
    return 0.0;
  }

  for (int i = 0; i < v1.getSize(); ++i)
    sum +=  v1[i] * v2[i];
  return sum;
}
