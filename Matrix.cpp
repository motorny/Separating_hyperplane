#include "Matrix.h"
#define EPS 1e-6
#include <iostream>
#include <iomanip>
using namespace std;
Matrix::Matrix()
{
}
Matrix::Matrix(const Matrix & m)
{
  row = m.getRowCnt();
  col = m.getColCnt();

  data = new mT*[row];
  for (int rI = 0; rI < row; ++rI)
  {
    data[rI] = new mT[col];
    for (int cI = 0; cI < col; ++cI)
      data[rI][cI] = m[rI][cI];
  }
}

Matrix::Matrix(int row_, int col_)
{
  row = row_;
  col = col_;
  data = new mT*[row];
  for (int rI = 0; rI < row; ++rI)
  {
    data[rI] = new mT[col];
    for (int cI = 0; cI < col; ++cI)
      data[rI][cI] = 0.0;
  }
}

Matrix::~Matrix()
{
  for (int rI = 0; rI < row; ++rI)
  {
    delete[] data[rI];
  }
  delete[] data;
}

Matrix & Matrix::operator=(const Matrix & m)
{
  if (row != m.getRowCnt() || col != m.getColCnt())
  {
    // if matrix has differtn size

    for (int rI = 0; rI < row; ++rI)
      delete[] data[rI];
    delete[] data;

    row = m.getRowCnt();
    col = m.getColCnt();

    data = new mT*[row];
    for (int rI = 0; rI < row; ++rI)
    {
      data[rI] = new mT[col];
      for (int cI = 0; cI < col; ++cI)
        data[rI][cI] = m[rI][cI];
    }
  }
  else
  {
    for (int rI = 0; rI < row; ++rI)
    {
      for (int cI = 0; cI < col; ++cI)
        data[rI][cI] = m[rI][cI];
    }
  }

  return (*this);
}

void Matrix::Extend(int addRowCnt, int addColCnt)
{
  int newRowCnt = row + addRowCnt;
  int newColCnt = col + addColCnt;

  mT** tempData = new mT*[newRowCnt];
  for (int rI = 0; rI < newRowCnt; ++rI)
  {
    tempData[rI] = new mT[newColCnt];
    if (rI >= row)
    {
      for (int cI = 0; cI < newColCnt; ++cI)
        tempData[rI][cI] = 0.0;
      continue;
    }
    // else - copy old values
    for (int cI = 0; cI < col; ++cI)
      tempData[rI][cI] = data[rI][cI];
    for (int cI = col; cI < newColCnt; ++cI)
      tempData[rI][cI] = 0.0;
  }

  for (int rI = 0; rI < row; ++rI)
    delete[] data[rI];
  delete[] data;

  data = tempData;
  row = newRowCnt;
  col = newColCnt;
}


Matrix  Matrix::MatrixMulMatrix(const IndexSet& iSetRowsA, const IndexSet& iSetColsB, const Matrix& mB) const
{
  int resRowCnt = iSetRowsA.GetSize();
  int resColCnt = iSetColsB.GetSize();

  Matrix res(resRowCnt, resColCnt);

  int resRowInd = 0, resColInd = 0;

  for (auto iRow : iSetRowsA.data)
  {
    resColInd = 0;
    for (auto iCol : iSetColsB.data)
    {
      mT sum = 0.0;
      for (int ind = 0; ind < col; ++ind)
      {
        sum += data[iRow][ind] * mB[ind][iCol];
      }
      res[resRowInd][resColInd] = sum;
      resColInd++;
    }
    resRowInd++;
  }

  return res;
}

void Matrix::Print(const char * s) const
{
  std::cout << s << std::endl;
  for (int iRow = 0; iRow < row; iRow++)
  {
    for (int iCol = 0; iCol < col; iCol++)
    {
      std::cout << std::setw(5) << data[iRow][iCol] << " ";
    }
    std::cout << "" << std::endl;
  }
}

void Matrix::Print(const char * s, const IndexSet & ColSet) const
{
  std::cout << s << std::endl;
  for (int iRow = 0; iRow < row; iRow++)
  {
    for (auto iCol : ColSet.data)
    {
      std::cout << std::setw(5) << data[iRow][iCol] << " ";
    }
    std::cout << "" << std::endl;
  }
}

bool Matrix::FullRang(const IndexSet & iSetRows, const IndexSet & iSetCols) const
{

  int rowCnt = iSetRows.GetSize();
  int colCnt = iSetCols.GetSize();

  Matrix m(rowCnt, colCnt);

  std::list<int>::const_iterator iterCol = iSetCols.data.cbegin();
  std::list<int>::const_iterator iterRow = iSetRows.data.cbegin();

  for (int iRow = 0; iRow < rowCnt; iRow++)
  {
    iterCol = iSetCols.data.cbegin();
    for (int iCol = 0; iCol < colCnt; iCol++)
    {
      m[iRow][iCol] = data[*iterRow][*iterCol];
      iterCol++;
    }
    iterRow++;
  }

  for (int iDiag = 0; iDiag < rowCnt - 1; iDiag++)
  {
    int maxInd = iDiag;
    double max = fabs(m[maxInd][iDiag]);


    for (int iRow = iDiag; iRow < rowCnt; iRow++)
    {
      if (max < fabs(m[iRow][iDiag]))
      {
        max = fabs(m[iRow][iDiag]);
        maxInd = iRow;
      }
    }

    if (max < EPS)
      return false;

    if (maxInd != iDiag)
    {
      for (int iColSw = 0; iColSw < colCnt; iColSw++)
      {
        double temp = m[maxInd][iColSw];
        m[maxInd][iColSw] = m[iDiag][iColSw];
        m[iDiag][iColSw] = temp;
      }
    }
    // m.Print("RowSwap:");


    for (int iRow = iDiag + 1; iRow < rowCnt; iRow++)
    {
      double multer = -m[iRow][iDiag] / m[iDiag][iDiag];
      for (int iColSummer = iDiag; iColSummer < colCnt; iColSummer++)
      {
        m[iRow][iColSummer] += m[iDiag][iColSummer] * multer;
      }
    }
    // m.Print("RowIter:");
  }

  //m.Print("Forward:");

  return fabs(m[rowCnt - 1][colCnt - 1]) > EPS;

}

mT * Matrix::operator[](size_t ind)
{
  return data[ind];
}

const mT * Matrix::operator[](size_t ind) const
{
  return data[ind];
}


Matrix Matrix::getInvertible(const IndexSet& iSetRows, const IndexSet& iSetCols) const
{
  int rowCnt = iSetRows.GetSize();
  int colCnt = iSetCols.GetSize() * 2;

  Matrix m(rowCnt, colCnt);

  std::list<int>::const_iterator iterCol = iSetCols.data.cbegin();
  std::list<int>::const_iterator iterRow = iSetRows.data.cbegin();

  for (int iRow = 0; iRow < rowCnt; iRow++)
  {
    iterCol = iSetCols.data.cbegin();
    for (int iCol = 0; iCol < colCnt / 2; iCol++)
    {
      m[iRow][iCol] = data[*iterRow][*iterCol];
      iterCol++;
    }
    iterRow++;
    m[iRow][colCnt / 2 + iRow] = 1.0;
  }

  //m.Print("Extended:");
  for (int iDiag = 0; iDiag < rowCnt - 1; iDiag++)
  {
    int maxInd = iDiag;
    double max = fabs(m[maxInd][iDiag]);


    for (int iRow = iDiag; iRow < rowCnt; iRow++)
    {
      if (max < fabs(m[iRow][iDiag]))
      {
        max = fabs(m[iRow][iDiag]);
        maxInd = iRow;
      }
    }

    if (maxInd != iDiag)
    {
      for (int iColSw = 0; iColSw < colCnt; iColSw++)
      {
        double temp = m[maxInd][iColSw];
        m[maxInd][iColSw] = m[iDiag][iColSw];
        m[iDiag][iColSw] = temp;
      }
    }
    //m.Print("RowSwap:");


    for (int iRow = iDiag + 1; iRow < rowCnt; iRow++)
    {
      double multer = -m[iRow][iDiag] / m[iDiag][iDiag];
      for (int iColSummer = iDiag; iColSummer < colCnt; iColSummer++)
      {
        m[iRow][iColSummer] += m[iDiag][iColSummer] * multer;
      }
    }
    //m.Print("RowIter:");
  }

  //m.Print("Forward:");
  for (int iDiag = rowCnt - 1; iDiag >= 0; iDiag--)
  {
    mT divider = m[iDiag][iDiag];
    for (int iColSummer = iDiag; iColSummer < colCnt; iColSummer++)
    {
      m[iDiag][iColSummer] /= divider;
    }
    // m.Print("Normalized:");
    for (int iRow = iDiag - 1; iRow >= 0; iRow--)
    {
      double multer = -m[iRow][iDiag];
      for (int iColSummer = iDiag; iColSummer < colCnt; iColSummer++)
      {
        m[iRow][iColSummer] += m[iDiag][iColSummer] * multer;
      }
    }
  }
  //m.Print("Backward:");

  Matrix res(rowCnt, colCnt / 2);
  for (int iRow = 0; iRow < rowCnt; iRow++)
  {
    for (int iCol = 0; iCol < colCnt / 2; iCol++) {
      res[iRow][iCol] = m[iRow][iCol + colCnt / 2];
    }
  }

  return res;
}

Matrix operator*(const Matrix & m1, const Matrix & m2)
{
  if (m1.getColCnt() != m2.getRowCnt())
  {
    std::cout << "MATRIX mul MATRIX DIMENSION ERROR" << std::endl;
    return Matrix();
  }

  Matrix res(m1.getRowCnt(), m2.getColCnt());

  for (int rI = 0; rI < m1.getRowCnt(); ++rI)
  {
    for (int cI = 0; cI < m2.getColCnt(); ++cI)
    {
      mT sum = 0.0;
      for (int ind = 0; ind < m1.getColCnt(); ++ind)
      {
        sum += m1[rI][ind] * m2[ind][cI];
      }
      res[rI][cI] = sum;
    }
  }

  return res;
}

Vec operator*(const Matrix & m, const Vec & v)
{
  if (m.getColCnt() != v.getSize())
  {
    cout << "MATRIX mul VEC DIMENSION ERROR" << endl;
    return Vec();
  }

  Vec res(m.getRowCnt());

  for (int rI = 0; rI < m.getRowCnt(); ++rI)
  {
    mT sum = 0.0;
    for (int ind = 0; ind < m.getColCnt(); ++ind)
    {
      sum += m[rI][ind] * v[ind];
    }
    res[rI] = sum;
  }

  return res;
}

Vec operator*(const Vec & v, const Matrix & m)
{
  if (v.getSize() != m.getRowCnt())
  {
    cout << "Vec mul Matrix DIMENSION ERROR" << endl;
    return Vec();
  }

  Vec res(m.getColCnt());

  for (int cI = 0; cI < m.getColCnt(); ++cI)
  {
    mT sum = 0.0;
    for (int ind = 0; ind < m.getRowCnt(); ++ind)
    {
      sum += m[ind][cI] * v[ind];
    }
    res[cI] = sum;
  }

  return res;
}
