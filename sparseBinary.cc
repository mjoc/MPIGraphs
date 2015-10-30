#include "SparseBinary.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <set>
#include <ctime>
#include "Runlog.hpp"


SparseBinary::SparseBinary(){
  _nRows = 0;
  _nCols = 0;
}

SparseBinary::SparseBinary(const SparseBinary& mat){
  _cols = mat._cols;
  _rowIndex = mat._rowIndex;
  _nRows = mat._nRows;
  _nCols = mat._nCols;
}


SparseBinary& SparseBinary::operator=(const SparseBinary& smat){
  SparseBinary temp(smat);

  std::swap(_nCols, temp._nCols);
  std::swap(_nRows, temp._nRows);
  std::swap(_cols, temp._cols);
  std::swap(_rowIndex, temp._rowIndex);
  return *this;
}


SparseBinary::SparseBinary(int nRows, int nCols, const std::vector<int>& colsArray, const std::vector<int>& rowStartIdxArray){
  _cols = colsArray;
  _rowIndex = rowStartIdxArray;
  _nRows = nRows;
  _nCols = nCols;
}


int SparseBinary::operator()(int iRow, int iCol){
  int result = 0;
  
  if(_rowIndex[iRow]!=_rowIndex[iRow+1]){
    for(int i = _rowIndex[iRow]; i < _rowIndex[iRow+1]; i++){
      if(_cols[i]==iCol){
        result = 1;
        break;
      }
    }
  }
  return result;
}

void SparseBinary::print(){
  int iRow, iCol, iCurrentCol, iData;
  
  if(_cols.size()==0){
    std::cout << "<<Empty " << _nRows << " x " << _nCols << " Matrix>>" << std::endl;
  }else{
    if(_nRows > 100 && _nCols > 100){
      std::cout << "Matrix is " << _nRows << " x " <<  _nCols << ", too big to print" << std::endl;
    }else{
      iData = 0;
      for(iRow = 0; iRow <_nRows; iRow++){
        std::cout << iRow << " : ";
        if(iData < (int)_cols.size() &&  _rowIndex[iRow] != _rowIndex[iRow+1]){
          iCurrentCol = _cols[iData];
          for(iCol = 0; iCol < _nCols; iCol++){
            if(iCol == iCurrentCol && iData < _rowIndex[iRow+1]){
              std::cout << "1 ";
              iData++;
              if(iData < (int)_cols.size()){
                iCurrentCol = _cols[iData];
              }
            }else{
              std::cout << "0 ";
            }
          }
        }else{
          for(iCol = 0; iCol <_nCols; iCol++){
            std::cout << " ";
          }
        }
        std::cout << std::endl;
      }
    }
  }
}

void SparseBinary::print(Runlog& runlog){
  int iRow, iCol, iCurrentCol, iData;
  std::string matrixLine;
  
  if(_cols.size()==0){
    matrixLine = std::string("<<Empty ") + runlog.to_string(_nRows) + std::string(" x ") + runlog.to_string(_nCols) + std::string(" Matrix>>");
    runlog.writelog(matrixLine);
  }else{
    if(_nRows > 100 && _nCols > 100){
      matrixLine = std::string("Matrix is") + runlog.to_string(_nRows) + std::string(" x ") + runlog.to_string(_nCols)+ std::string(", too big to print");
      runlog.writelog(matrixLine);
    }else{
      iData = 0;
      if(_nRows> 0){
				for(iRow = 0; iRow <_nRows; iRow++){
          matrixLine = runlog.to_string(iRow) + std::string(" : ");

					if(iData < (int)_cols.size() &&  _rowIndex[iRow] != _rowIndex[iRow+1]){
						iCurrentCol = _cols[iData];

						for(iCol = 0; iCol < _nCols; iCol++){
							if(iCol == iCurrentCol && iData < _rowIndex[iRow+1]){
                matrixLine += std::string("1 ");

								iData++;
								if(iData < (int)_cols.size()){
									iCurrentCol = _cols[iData];
								}
							}else{
                matrixLine += std::string("0 ");

							}
						}
					}else{
						for(iCol = 0; iCol <_nCols; iCol++){
              matrixLine += std::string("0 ");

						}
					}
          runlog.writelog(matrixLine);

				}
      }
    }
  }
  return;
}
