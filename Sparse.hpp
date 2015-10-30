#ifndef __MOC_Sparse_
#define __MOC_Sparse_

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <set>
#include <ctime>
#include <stdexcept>
#include "Runlog.hpp"

template<typename DataType>
class Sparse;

template<typename DataType>
Sparse<DataType> operator-(const Sparse<DataType> &A, const Sparse<DataType> &B);

template<typename DataType>
std::vector<double> operator*(const Sparse<DataType> &A, const std::vector<double> &b);

template<typename DataType>
class Sparse{
public:
  int _nRows;
  int _nCols;
  std::vector<DataType> _data;
  std::vector<int> _cols;
  std::vector<int> _rowIndex;

  
  ~Sparse();
  Sparse();
  Sparse(int nRows, int nCols, const std::vector<DataType> &dataArray, const std::vector<int>& colsArray, const std::vector<int>& rowStartIdxArray);
  Sparse(const std::vector<DataType> &diag, const std::vector<DataType>& offDiag);
  Sparse(const Sparse&);
  Sparse& operator=(const Sparse&);
  void print() const;
  void print(Runlog& runlog) const;
  bool readFromFile(char *filename);
  int nRows()const {return _nRows;};
  int nCols()const {return _nCols;};
  int nData()const {return _rowIndex[_nRows];}
  DataType operator()(int iRow, int iCol);
  
  template<typename DataType2>
  friend std::vector<double> postMultBy(const Sparse<DataType2>& mat, const std::vector<double>& b);
};

template<typename DataType>
Sparse<DataType>::~Sparse(){

}

template<typename DataType>
Sparse<DataType>::Sparse(){
  _nRows = 0;
  _nCols = 0;
}


template<typename DataType>
Sparse<DataType>::Sparse(const Sparse& mat){
  _data = mat._data;
  _cols = mat._cols;
  _rowIndex = mat._rowIndex;
  _nRows = mat._nRows;
  _nCols = mat._nCols;
}

template<typename DataType>
Sparse<DataType>& Sparse<DataType>::operator=(const Sparse& smat){
  Sparse temp(smat);

  std::swap(_nCols, temp._nCols);
  std::swap(_nRows, temp._nRows);
  std::swap(_data, temp._data);
  std::swap(_cols, temp._cols);
  std::swap(_rowIndex, temp._rowIndex);
  return *this;
  
}

/*
template<typename DataType>
Sparse<DataType>::Sparse( int nRows, int nCols, int nData, DataType *dataArray, int *colsArray, int *rowStartIdxArray){
  bool success =true;
  //cout << "Number of data are " << nData << endl;
  try{
    _data = new DataType[nData];
    _cols = new int[nData];
    _rowIndex = new int[nRows+1];
  }catch(...){
    success = false;
  }
  if(success){
    _nData = nData;
    _nRows = nRows;
    _nCols = nCols;
    
    for(int i = 0; i < nData; i++){
      _data[i] =  dataArray[i];
      _cols[i] = colsArray[i];
    }
    for(int i = 0; i <= nRows; i++){
      _rowIndex[i] = rowStartIdxArray[i];
    }
    
  }
}
*/

template<typename DataType>
Sparse<DataType>::Sparse(int nRows, int nCols, const std::vector<DataType> &dataArray, const std::vector<int>& colsArray, const std::vector<int>& rowStartIdxArray){

  //cout << "Number of data are " << nData << endl;
  _data = dataArray;
  _cols = colsArray;
  _rowIndex = rowStartIdxArray;
  _nRows = nRows;
  _nCols = nCols;
  
}

template<typename DataType>
Sparse<DataType>::Sparse(const std::vector<DataType> &diag, const std::vector<DataType>& offDiag){
  _data.resize(diag.size()+2*offDiag.size());
  _cols.resize(diag.size()+2*offDiag.size());
  _rowIndex.resize(diag.size()+1);
  _nRows = (int)diag.size();
  _nCols = (int)diag.size();
  
  int i= 0;
  _rowIndex[i] = 0;
  _data[i*3] = diag[i];
  _data[i*3+1] = offDiag[i];
  _cols[i*3] = i;
  _cols[i*3+1] = i+1;
  
  for(i = 1; i < (int)diag.size()-1;i++){
    _rowIndex[i] = i*3 -1;
    _data[i*3-1] = offDiag[i-1];
    _data[i*3] = diag[i];
    _data[i*3+1] = offDiag[i];
    _cols[i*3-1] = i-1;
    _cols[i*3] = i;
    _cols[i*3+1] = i+1;
    
  }
  _rowIndex[i] = i*3 -1;
  _data[i*3-1] = offDiag[i-1];
  _data[i*3] = diag[i];
  _cols[i*3-1] = i-1;
  _cols[i*3] = i;
  i++;
  _rowIndex[i] =  i*3 -2;
  
}

template<typename DataType>
bool Sparse<DataType>::readFromFile(char *filename){
  bool allOk = true;
  
  _data.resize(0);
  _cols;
  _rowIndex;
  _nRows = 0;
  _nCols = 0;
  
  std::ifstream infile(filename,std::ios_base::in);
  //cout << "Reading in data from file " << filename << endl;
  std::string line = "", numberString;
  std::getline(infile, line);
  
  int nData, nRows, nCols, n1, n2;
  n1 = (int)line.find(' ');
  n2 = (int)line.find(' ',n1+1);
  numberString = line.substr (0,n1);
  nRows = atoi(numberString.c_str());
  numberString = line.substr (n1+1,n2-n1+1);
  nCols = atoi(numberString.c_str());
  //cout << "Matrix is " <<  nRows << " x " << nCols << endl;
  numberString = line.substr (n2+1);
  nData = atoi(numberString.c_str());
  //cout << "There are " << nData << " data entries" << endl;
  
  int *data = new int[nData];
  int *colNumber = new int[nData];
  int *rowIndex = new int[nRows+1];
  int i = 0, iRow, iCurrentRow = 0;
  rowIndex[0] = 0;
  while (getline(infile, line)){
    //break on a blank line
    if(line.find_first_not_of(' ') == std::string::npos){
      break;
    }else{
      n1 = line.find(' ');
      n2 = line.find(' ',n1+1);
      numberString = line.substr (0,n1);
      iRow = atoi(numberString.c_str());
      
      if(iRow - iCurrentRow > 0){
        if(iRow - iCurrentRow == 1){
          iCurrentRow++;
        }else{
          iCurrentRow++;
          while(iCurrentRow < iRow){
            //cout << "currentRow++ to " << iCurrentRow << endl;
            //cout << "+Adding " << i << " to " << iCurrentRow << " place in the array" << endl;
            rowIndex[iCurrentRow] = i;
            iCurrentRow++;
          }
        }
        //cout << " Adding " << i << " to " << iRow << " place in the array" << endl;
        rowIndex[iRow] = i;
      }
      numberString = line.substr (n1+1,n2-n1+1);
      colNumber[i] = atoi(numberString.c_str());
      numberString = line.substr (n2+1);
      data[i] = atoi(numberString.c_str());
      i++;
    }
  }
  while(iRow < nRows){
    iRow++;
    rowIndex[iRow] = nData;
    
  }
  
  if(allOk){
    _nRows = nRows;
    _nCols = nCols;
    _data = data;
    _cols = colNumber;
    _rowIndex = rowIndex;
    
  }else{
    _nRows = 0;
    _nCols = 0;
    _data.resize(0);
    _cols.resize(0);
    _rowIndex.resize(0);
  }
  
  
  return allOk;
}

template<typename DataType>
DataType Sparse<DataType>::operator()(int iRow, int iCol){
  DataType result = (DataType)0;
  
  if(_rowIndex[iRow]!=_rowIndex[iRow+1]){
    for(int i = _rowIndex[iRow]; i < _rowIndex[iRow+1]; i++){
      if(_cols[i]==iCol){
        result = _data[i];
        break;
      }
    }
  }
  return result;
}

template<typename DataType>
void Sparse<DataType>::print() const{
  int iRow, iCol, iCurrentCol, iData;
  /*
  std::cout << "Data: ";
  for(size_t i= 0; i < _data.size();i++){
    std::cout << _data[i] << " ";
  }
  std::cout << std::endl;
  
  std::cout << "Cols: ";
  for(size_t i= 0; i < _cols.size();i++){
    std::cout << _cols[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "RowIndex: ";
  for(size_t i= 0; i < _rowIndex.size();i++){
    std::cout << _rowIndex[i] << " ";
  }
  std::cout << std::endl;
*/
  if(_data.size()==0){
    std::cout << "<<Empty " << _nRows << " x " << _nCols << " Matrix>>" << std::endl;
  }else{
    if(_nRows > 100 && _nCols > 100){
      std::cout << "Matrix is " << _nRows << " x " <<  _nCols << ", too big to print" << std::endl;
    }else{
      //cout << "Printing Sparse Mat" << endl;
      iData = 0;
      for(iRow = 0; iRow <_nRows; iRow++){
        std::cout << iRow << " : ";
        if(iData < (int)_data.size() &&  _rowIndex[iRow] != _rowIndex[iRow+1]){
          iCurrentCol = _cols[iData];
          //cout << "there are data here" << endl;
          for(iCol = 0; iCol < _nCols; iCol++){
            if(iCol == iCurrentCol && iData < _rowIndex[iRow+1]){
              std::cout << _data[iData] << " ";
              iData++;
              if(iData < (int)_data.size()){
                iCurrentCol = _cols[iData];
              }
            }else{
              std::cout << 0 << " ";
            }
          }
        }else{
          for(iCol = 0; iCol <_nCols; iCol++){
            std::cout << 0 << " ";
          }
        }
        std::cout << std::endl;
      }
    }
  }
}

template<typename DataType>
void Sparse<DataType>::print(Runlog& runlog) const{
  int iRow, iCol, iCurrentCol, iData;
  std::string matrixLine;
  
  if(_data.size()==0){
    matrixLine = std::string("<<Empty ") + runlog.to_string(_nRows) + std::string(" x ") + runlog.to_string(_nCols) + std::string(" Matrix>>");
    runlog.writelog(matrixLine);
    //std::cout << "<<Empty " << _nRows << " x " << _nCols << " Matrix>>" << std::endl;
  }else{
    if(_nRows > 100 && _nCols > 100){
      matrixLine = std::string("Matrix is") + runlog.to_string(_nRows) + std::string(" x ") + runlog.to_string(_nCols)+ std::string(", too big to print");
      runlog.writelog(matrixLine);
      //std::cout << "Matrix is " << _nRows << " x " <<  _nCols << ", too big to print" << std::endl;
    }else{
      //cout << "Printing Sparse Mat" << endl;
      iData = 0;
      if(_nRows> 0){
				for(iRow = 0; iRow <_nRows; iRow++){
          matrixLine = runlog.to_string(iRow) + std::string(" : ");
          //std::cout << iRow << " : ";
					if(iData < (int)_data.size() &&  _rowIndex[iRow] != _rowIndex[iRow+1]){
						iCurrentCol = _cols[iData];
						//cout << "there are data here" << endl;
						for(iCol = 0; iCol < _nCols; iCol++){
							if(iCol == iCurrentCol && iData < _rowIndex[iRow+1]){
                matrixLine += runlog.to_string(_data[iData]) + std::string(" ");
                //std::cout << _data[iData] << " ";
								iData++;
								if(iData < (int)_data.size()){
									iCurrentCol = _cols[iData];
								}
							}else{
                matrixLine += std::string("0 ");
                //std::cout << 0 << " ";
							}
						}
					}else{
						for(iCol = 0; iCol <_nCols; iCol++){
              matrixLine += std::string("0 ");
              //std::cout << 0 << " ";
						}
					}
          runlog.writelog(matrixLine);
          //std::cout << std::endl;
				}
      }
    }
  }
  return;
}


template<typename DataType>
Sparse<DataType> operator-(const Sparse<DataType> &Amat, const Sparse<DataType> &Bmat)
{
  if(Amat._nRows != Bmat._nRows || Amat._nCols != Bmat._nCols){
    throw std::logic_error("Rows and Columns of Sparse matrix substraction are not consistent");
  }
  
  std::vector<DataType> data;
  std::vector<int> cols;
  std::vector<int> rowIndex(Amat._nRows+1);
  std::vector< std::set<int> > combinedDataPositions(Amat._nRows);
  
  
  std::clock_t start = std::clock();
  
  
  
  int nData = 0;
  std::set<int> columns;
  for(int i = 0; i < Amat._nRows; i++){
    for(int j = Amat._rowIndex[i]; j < Amat._rowIndex[i+1];j++){
      combinedDataPositions[i].insert(Amat._cols[j]);
    }
    for(int j = Bmat._rowIndex[i]; j < Bmat._rowIndex[i+1];j++){
      combinedDataPositions[i].insert(Bmat._cols[j]);
    }
    nData += (int)combinedDataPositions[i].size();
  }
  
  //std::cout << "Combined data positions in " << (std::clock() - start) / (double)(CLOCKS_PER_SEC) << " seconds" << std::endl;
  
  std::stringstream message;
  
  
  start = std::clock();
  
  std::set<int>::iterator it;
  data.resize(nData);
  cols.resize(nData);
  int  iA = 0, iB = 0, iAB = 0;
  for(int iRow = 0; iRow < Amat._nRows; iRow++){
    for (it=combinedDataPositions[iRow].begin(); it!=combinedDataPositions[iRow].end(); ++it){
      if(iA < Amat._rowIndex[iRow+1] && Amat._cols[iA] == *it){
        data[iAB] += Amat._data[iA];
        cols[iAB] += Amat._cols[iA];
        iA++;
      }
      
      if(iB < Bmat._rowIndex[iRow+1] && Bmat._cols[iB] == *it){
        data[iAB] -= Bmat._data[iB];
        cols[iAB] += Bmat._cols[iB];
        iB++;
      }
      iAB++;
    }
    rowIndex[iRow+1] = iAB;
  }
  
  std::cout << "Subtraction in " << (std::clock() - start) / (double)(CLOCKS_PER_SEC) << " seconds" << std::endl;
  
  Sparse<DataType> result( Amat._nRows, Amat._nCols, data, cols, rowIndex);
  
  return result;
  
}





//template<typename DataType>
//Sparse<DataType> operator-(const Sparse<DataType> &Amat, const Sparse<DataType> &Bmat)
//{
//  if(Amat._nRows != Bmat._nRows || Amat._nCols != Bmat._nCols){
//      throw MyException("Rows and Columns of Sparse matrix substraction are not consistent");
//  }
//  
//  std::vector<DataType> data;
//  std::vector<int> cols;
//  std::vector<int> rowIndex(Amat._nRows+1);
//  std::vector< std::set<int> > combinedDataPositions(Amat._nRows);
//  
//  
//  std::clock_t start = std::clock();
//    
//  
//  
//  int nData = 0;
//  std::set<int> columns;
//  for(int i = 0; i < Amat._nRows; i++){
//    for(int j = Amat._rowIndex[i]; j < Amat._rowIndex[i+1];j++){
//      combinedDataPositions[i].insert(Amat._cols[j]);
//    }
//    for(int j = Bmat._rowIndex[i]; j < Bmat._rowIndex[i+1];j++){
//      combinedDataPositions[i].insert(Bmat._cols[j]);
//    }
//    nData += (int)combinedDataPositions[i].size();
//  }
//  
//  std::cout << "Combined data positions in " << (std::clock() - start) / (double)(CLOCKS_PER_SEC) << " seconds" << std::endl;
//  
//  std::stringstream message;
//  
//  
//  start = std::clock();
//  
//  std::set<int>::iterator it;
//  data.resize(nData);
//  cols.resize(nData);
//  int  iA = 0, iB = 0, iAB = 0;
//  for(int iRow = 0; iRow < Amat._nRows; iRow++){
//    for (it=combinedDataPositions[iRow].begin(); it!=combinedDataPositions[iRow].end(); ++it){
//      if(iA < Amat._rowIndex[iRow+1] && Amat._cols[iA] == *it){
//        data[iAB] += Amat._data[iA];
//        cols[iAB] += Amat._cols[iA];
//        iA++;
//      }
//
//      if(iB < Bmat._rowIndex[iRow+1] && Bmat._cols[iB] == *it){
//        data[iAB] -= Bmat._data[iB];
//        cols[iAB] += Bmat._cols[iB];
//        iB++;
//      }
//      iAB++;
//    }
//    rowIndex[iRow+1] = iAB;
//  }
//  
//  std::cout << "Subtraction in " << (std::clock() - start) / (double)(CLOCKS_PER_SEC) << " seconds" << std::endl;
//  
//  Sparse<DataType> result( Amat._nRows, Amat._nCols, data, cols, rowIndex);
//  
//  return result;
//  
//}

//template<typename DataType>
//Sparse<DataType> operator-(const Sparse<DataType> &Amat, const Sparse<DataType> &Bmat)
//{
//  if(Amat._nRows != Bmat._nRows || Amat._nCols != Bmat._nCols){
//    throw MyException("Rows and Columns of Sparse matrix substraction are not consistent");
//  }
//  
//  std::cout << "¢¢¢¢¢¢¢¢¢¢here" << std::endl;
//  
//  int iData = -1, jData = -1, iIndex = 0, jIndex = -1;
//  int iRow = 0, jRow = 0, nCols=Amat._nCols;
//  //int commonCells = 0;
//  bool notFinishedA = true, notFinishedB = true;
//  bool waitB = false;
//
//	if(Amat._nData == 0){
//		notFinishedA = false;
//	}
//	if(Bmat._nData == 0){
//		notFinishedB = false;
//	}
//  
//  //int newData[maxData];
//  std::vector<int> newDataIndex(Amat._nData +  Bmat._nData);
//  std::vector<int> newDataOrigin(Amat._nData +  Bmat._nData);
//  DataType newData[Amat._nData +  Bmat._nData];
//  std::vector<int> newCols(Amat._nData +  Bmat._nData);
//  int iSmallData=0, iSmallLastDouble = -1, nSinglesToPutIn = 0, iSmallRow = 0;
//  std::vector<DataType> smallData(Amat._nData +  Bmat._nData);
//  std::vector<int> smallDataIndex(Amat._nData +  Bmat._nData);
//  std::vector<int> smallCols(Amat._nData +  Bmat._nData);
//  std::vector<int> smallRows(Amat._nData +  Bmat._nData);
//  std::vector<int> smallRowIndex(Amat._nRows+1);
//  smallRowIndex[0] = 0;
//  int iNewData = 0;
//  int iSearch, jSearch;
//
//  
//  while(notFinishedA || notFinishedB){
//    if(notFinishedA){
//      iData++;
//      
//      while(iData >= Amat._rowIndex[iRow+1]){
//        iRow++;
//      }
//      
//      //while(iData == Amat._rowIndex[iRow] && iData == Amat._rowIndex[iRow+1]){
//      //  iRow++;
//      //}
//      
//      iIndex = iRow * Amat._nCols + Amat._cols[iData];
//      
//      iSearch = iNewData;
//      while(iSearch > 0 &&  newDataIndex[iSearch-1] > iIndex){
//        iSearch--;
//      }
//      
//      if(iSearch < iNewData){
//        for(int i = iNewData; i > iSearch; i--){
//          newDataIndex[i] = newDataIndex[i-1];
//          newDataOrigin[i] = newDataOrigin[i-1];
//          newData[i] = newData[i-1];
//          newCols[i] = newCols[i-1];
//        }
//      }
//      newDataIndex[iSearch] = iIndex;
//      newDataOrigin[iSearch] = 1;
//      newData[iSearch] = Amat._data[iData];
//      newCols[iSearch] = Amat._cols[iData];
//      
//      // When we reach data in both matrices we then fill in teh fully compressed
//      // data to that point
//      if(iSearch > 0 && (newDataIndex[iSearch]==newDataIndex[iSearch-1])){
//        nSinglesToPutIn = iNewData - iSmallLastDouble -2;
//        if(nSinglesToPutIn > 0){
//          //cout << "A There are " << nSinglesToPutIn  << " singles" << endl;
//          for(int i = iSmallLastDouble+1; i < iSmallLastDouble+nSinglesToPutIn+1;i++){
//            smallData[iSmallData] = newData[i];
//            smallCols[iSmallData] = newCols[i];
//            smallRows[iSmallData] = newDataIndex[i]/nCols;
//            if(iSmallData > 0 &&  smallRows[iSmallData] != smallRows[iSmallData-1]){
//              //cout << "New row " << smallRows[iSmallData-1] << " to " << smallRows[iSmallData] << endl;
//              while(iSmallRow < smallRows[iSmallData]){
//                iSmallRow++;
//                smallRowIndex[iSmallRow] = iSmallData;
//              }
//            }
//            
//            smallDataIndex[iSmallData] = newDataIndex[i];
//            iSmallData++;
//          }
//        }
//        smallData[iSmallData] = newData[iSearch-1] + newData[iSearch];
//        smallCols[iSmallData] = newCols[iSearch];
//        smallRows[iSmallData] = newDataIndex[iSearch]/nCols;
//        if(iSmallData > 0 &&  smallRows[iSmallData] != smallRows[iSmallData-1]){
//          //cout << "New row " << smallRows[iSmallData-1] << " to " << smallRows[iSmallData] << endl;
//          while(iSmallRow < smallRows[iSmallData]){
//            iSmallRow++;
//            smallRowIndex[iSmallRow] = iSmallData;
//          }
//        }
//        smallDataIndex[iSmallData] = newDataIndex[iSearch];
//        iSmallLastDouble = iNewData;
//        iSmallData++;
//      }
//      
//      iNewData++;
//      
//      if(iData == Amat._nData-1){
//        notFinishedA = false;
//      }
//    }
//    if(iIndex > jIndex){
//      waitB = false;
//    }
//    
//    while(notFinishedB && !waitB){
//      jData++;
//      
//      while(jData >= Bmat._rowIndex[jRow+1]){
//        jRow++;
//      }
//      while(jData == Bmat._rowIndex[jRow] && jData == Bmat._rowIndex[jRow+1]){
//        jRow++;
//      }
//      jIndex = jRow * Bmat._nCols + Bmat._cols[jData];
//      
//      jSearch = iNewData;
//      while(jSearch > 0 &&  newDataIndex[jSearch-1] > jIndex){
//        jSearch--;
//      }
//      if(jSearch < iNewData){
//        for(int i = iNewData; i > jSearch; i--){
//          newDataIndex[i] = newDataIndex[i-1];
//          newDataOrigin[i] = newDataOrigin[i-1];
//          newData[i] = newData[i-1];
//          newCols[i] = newCols[i-1];
//        }
//      }
//      newDataIndex[jSearch] = jIndex;
//      newDataOrigin[jSearch] = 2;
//      newData[jSearch] = -Bmat._data[jData];
//      newCols[jSearch] = Bmat._cols[jData];
//      
//      // When we reach data in both matrices we then fill in teh fully compressed
//      // data to that point
//      if(jSearch > 0 && newDataIndex[jSearch] == newDataIndex[jSearch-1]){
//        nSinglesToPutIn = iNewData - iSmallLastDouble -2;
//        if(nSinglesToPutIn > 0){
//          //cout << "B There are " << nSinglesToPutIn  << " singles" << endl;
//          for(int i = iSmallLastDouble+1; i < iSmallLastDouble+nSinglesToPutIn+1;i++){
//            smallData[iSmallData] = newData[i];
//            smallCols[iSmallData] = newCols[i];
//            smallRows[iSmallData] = newDataIndex[i]/nCols;
//            if(iSmallData > 0 &&  smallRows[iSmallData] != smallRows[iSmallData-1]){
//              //cout << "New row " << smallRows[iSmallData-1] << " to " << smallRows[iSmallData] << endl;
//              while(iSmallRow < smallRows[iSmallData]){
//                iSmallRow++;
//                smallRowIndex[iSmallRow] = iSmallData;
//                
//              }
//            }
//            smallDataIndex[iSmallData] = newDataIndex[i];
//            iSmallData++;
//          }
//        }
//        smallData[iSmallData] = newData[jSearch-1] + newData[jSearch];
//        smallCols[iSmallData] = newCols[jSearch-1];
//        smallRows[iSmallData] = newDataIndex[jSearch]/nCols;
//        if(iSmallData > 0 &&  smallRows[iSmallData] != smallRows[iSmallData-1]){
//          //cout << "New row " << smallRows[iSmallData-1] << " to " << smallRows[iSmallData] << endl;
//          while(iSmallRow < smallRows[iSmallData]){
//            iSmallRow++;
//            smallRowIndex[iSmallRow] = iSmallData;
//            
//          }
//        }
//        smallDataIndex[iSmallData] = newDataIndex[jSearch];
//        iSmallLastDouble = iNewData;
//        iSmallData++;
//      }
//      iNewData++;
//      
//      if(jIndex >= iIndex){
//        waitB = true;
//      }else{
//        waitB = false;
//      }
//      
//      if(jData == Bmat._nData-1){
//        notFinishedB = false;
//      }
//    }
//    
//  }
//  
//  nSinglesToPutIn = iNewData -1 - iSmallLastDouble ; // (iNewData -1 is the last data point)
//  if(nSinglesToPutIn > 0){
//    //cout << "C There are " << nSinglesToPutIn  << " singles" << endl;
//    for(int i = iSmallLastDouble+1; i < iSmallLastDouble+nSinglesToPutIn+1;i++){
//      smallData[iSmallData] = newData[i];
//      smallCols[iSmallData] = newCols[i];
//      smallRows[iSmallData] = newDataIndex[i]/nCols;
//      if(iSmallData > 0 &&  smallRows[iSmallData] != smallRows[iSmallData-1]){
//        //cout << "New row " << smallRows[iSmallData-1] << " to " << smallRows[iSmallData] << endl;
//        while(iSmallRow < smallRows[iSmallData]){
//          iSmallRow++;
//          smallRowIndex[iSmallRow] = iSmallData;
//          
//        }
//      }
//      smallDataIndex[iSmallData] = newDataIndex[i];
//      iSmallData++;
//    }
//  }
//  smallRowIndex[Amat._nRows] = iSmallData;
//  //cout << "There are " << commonCells << " common cells" << endl;
//  
//  //
//  //for(int i = 0; i < iNewData; i++){
//  //  cout << newDataIndex[i] << " " << newDataOrigin[i]  << " " << newData[i] << endl;
//  //}
//  //cout << "Compressed" << endl;
//  //for(int i = 0; i < iSmallData; i++){
//  //  cout << smallDataIndex[i] << " " << smallData[i] << " " << smallRows[i] << " " << smallCols[i] << endl;
//  //}
//  //cout << "Row pointer:" << endl;
//  //for(int i = 0; i <= Amat._nCols; i++){
//  //  cout << smallRowIndex[i] << " " ;
//  //}
//  //cout << endl;
//  //int data[A._nData + B._nData - commonCells];
//  //int cols[A._nData + B._nData - commonCells];
//  //int rowIndex[A._nData + B._nData - commonCells];
//  
//  std::cout << "$$$$$$$$There" << std::endl;
//  
//  Sparse<DataType> result( Amat._nRows, Amat._nCols, iNewData, smallData, smallCols, smallRowIndex);
//  
//  return result;
//}



template<typename DataType>
std::vector<double> postMultBy(const Sparse<DataType>& mat, const std::vector<double>& b){
  if(b.size() != mat._nCols){
    throw std::logic_error("Rows and Columns of Sparse matrix / My_vector multiplication are not consistent");
  }
  
  std::vector<double> result(mat._nCols);
  
  double sumproduct;
  int iDataIndex, iCol;
  for(int iRow = 0; iRow < b.size(); iRow++){
    if(mat._rowIndex[iRow] < mat._rowIndex[iRow+1]){
      sumproduct = 0.0;
      for(iDataIndex = mat._rowIndex[iRow]; iDataIndex< mat._rowIndex[iRow+1];  iDataIndex++){
        iCol = mat._cols[iDataIndex];
        sumproduct += (double)mat._data[iDataIndex] * b[iCol];
      }
      result[iRow] = sumproduct;
    }
  }
  
  return result;
}

#endif
