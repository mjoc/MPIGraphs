#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <ctime>
#include <stdexcept>
#include "Graph.hpp"
#include "Sparse.hpp"
#include "SparseBinary.hpp"


bool AdjSMatGraph::readFromFile(const char *filename, int formatCode){
  bool allOk = true;
  int nVertices = 0, nEdges = 0;
  std::ifstream infile(filename);

  std::cout << "Reading in data from file " << filename << std::endl;
  std::string line = "";
  allOk = false;
  bool oneBased = false;
  while (getline(infile, line)){
    //break on a blank line
    if(line.find_first_not_of(' ') == std::string::npos){
      break;
    }else{
      allOk = true;
      nVertices++;
      std::stringstream strstr(line);
      std::string word = "";
      switch(formatCode){
        case 0:
          getline(strstr,word, ';');
          if(nVertices == 1){
            if(atoi(word.c_str()) == 1){
              oneBased = true;
            }
          }
          while (getline(strstr,word, ';')){
            nEdges++;
            //std::cout << " " << iCol << "(" << iRow*nVertices+iCol << ")" << " ";
          }
          break;
        case 1:
          getline(strstr,word, ':');
          if(nVertices == 1){
            if(atoi(word.c_str()) == 1){
              oneBased = true;
            }
          }
          getline(strstr,word, ' ');
          while (getline(strstr,word, ' ')){
            nEdges++;
            //std::cout << " " << iCol << "(" << iRow*nVertices+iCol << ")" << " ";
          }
          break;
        default:
          throw std::logic_error("Format code for graph file not found\n");
      }
    }
  }
  if(oneBased){
    std::cout << "One (1) based index, changing to zero based\n";
  };
  std::cout << "Number of vertices: " << nVertices << std::endl;
  std::cout << "Number of edges: " << nEdges << std::endl;
  if(allOk){
    _nRows = nVertices;
    _nCols = nVertices;
    
    _cols.resize(nEdges);
    _rowIndex.resize(nVertices+1);
    
    infile.clear(); //clear error flags
    infile.seekg(0, std::ios::beg); //set the file get pointer back to the beginning
    
    int iRow, iCol, iLine, iCurrentRow = 0, iEdge = 0;
    std::vector<int> colsInRow;

    for(iLine = 0 ; iLine < nVertices; iLine++){
      getline(infile, line);
      //break on a blank line
      if(line.find_first_not_of(' ') == std::string::npos)
      {
        break;
      }
      std::stringstream strstr(line);
      std::string word = "";
      switch(formatCode){
        case 0:
          getline(strstr,word, ';');
          iRow = atoi(word.c_str());
          if(oneBased){iRow--;}
          
          while(iCurrentRow < iRow){
            _rowIndex[iCurrentRow]= iEdge;
            iCurrentRow++;
          }

          _rowIndex[iCurrentRow]= iEdge;

          colsInRow.resize(0);
          while (getline(strstr,word, ';')){
            iCol = atoi(word.c_str());
            if(oneBased){iCol--;}
            colsInRow.push_back(iCol);
          }
          break;
        case 1:
          getline(strstr,word, ':');
          iRow = atoi(word.c_str());
          if(oneBased){iRow--;}
          while(iCurrentRow < iRow){
            _rowIndex[iCurrentRow]= iEdge;
            iCurrentRow++;
          }
          _rowIndex[iCurrentRow]= iEdge;
          getline(strstr,word, ' ');
          colsInRow.resize(0);
          while (getline(strstr,word, ' ')){
            iCol = atoi(word.c_str());
            if(oneBased){iCol--;}
            colsInRow.push_back(iCol);
          }
          break;
          
      }
      std::sort(colsInRow.begin(),colsInRow.end());
      for(std::vector<int>::iterator it = colsInRow.begin(); it != colsInRow.end(); ++it) {
        _cols[iEdge] = *it;
        iEdge++;
      }
      iCurrentRow++;
    }
    _rowIndex[nVertices] =  nEdges;

  }
  infile.close();
  return allOk;
}



bool AdjSMatGraph::readSloppyFromFile(char *filename, int formatCode){
  bool allOk = true;
  std::ifstream infile(filename);
  
  std::cout << "Reading in data from file " << filename << std::endl;
  std::string line = "";
  allOk = false;
  int maxVertex = 0;
  
  std::vector<std::set<int> > tempData;
  
  while (getline(infile, line)){
    //break on a blank line
    if(line.find_first_not_of(' ') == std::string::npos){
      break;
    }else{
      allOk = true;
      int iRow, iCol;
      std::stringstream strstr(line);
      std::string word = "";
      switch(formatCode){
        case 0:
          getline(strstr,word, ';');
          iRow = atoi(word.c_str());
          if(iRow > maxVertex){maxVertex = iRow;}
          while (getline(strstr,word, ';')){
            iCol = atoi(word.c_str());
            if(iCol > maxVertex){maxVertex = iCol;}
          }
          break;
        case 1:
          getline(strstr,word, ':');
          getline(strstr,word, ' ');
          while (getline(strstr,word, ' ')){
            iCol = atoi(word.c_str());
            if(iCol > maxVertex){maxVertex = iCol;}
          }
          break;
        default:
          throw std::logic_error("Format code for graph file not found\n");
      }
    }
  }
  std::cout << "Max vertex (0 based index) is " << maxVertex << std::endl;
  int nVertices = maxVertex + 1;
  
  if(allOk){
    tempData.resize(maxVertex+1);
 
    infile.clear(); //clear error flags
    infile.seekg(0, std::ios::beg); //set the file get pointer back to the beginning
    
    int iRow, iCol;
    std::vector<int> colsInRow;
    //std::cout << "check ";
    while (getline(infile, line)){
      //break on a blank line
      if(line.find_first_not_of(' ') == std::string::npos)
      {
        break;
      }
      std::stringstream strstr(line);
      std::string word = "";
      switch(formatCode){
        case 0:
          getline(strstr,word, ';');
          iRow = atoi(word.c_str());
          
          colsInRow.resize(0);
          while (getline(strstr,word, ';')){
            iCol = atoi(word.c_str());
            colsInRow.push_back(iCol);
            tempData[iRow].insert(iCol);
            tempData[iCol].insert(iRow);
          }
          break;
        case 1:
          getline(strstr,word, ':');
          iRow = atoi(word.c_str())-1;
          
          getline(strstr,word, ' ');
          colsInRow.resize(0);
          while (getline(strstr,word, ' ')){
            iCol = atoi(word.c_str())-1;
            colsInRow.push_back(iCol);
            tempData[iRow].insert(iCol);
            tempData[iCol].insert(iRow);
          }
          break;
          
      }
    }
    
    int nAllEdges = 0;
    for(std::vector<std::set<int> >::const_iterator itVertex = tempData.begin(); itVertex != tempData.end(); itVertex++){
      nAllEdges += (int)(*itVertex).size();
    }
    std::cout << "There are " << nAllEdges << " edges\n";
    
    _cols.resize(nAllEdges);
    _rowIndex.resize(nVertices+1);
    
    int iVertex = 0, iData = 0;
    _rowIndex[0] = 0;
    for(std::vector<std::set<int> >::const_iterator itVertex = tempData.begin(); itVertex != tempData.end(); itVertex++){
      for (std::set<int>::iterator itEdgeSet = (*itVertex).begin(); itEdgeSet != (*itVertex).end(); itEdgeSet++) {
        _cols[iData] = *itEdgeSet;
        iData++;
      }
      _rowIndex[iVertex+1] = iData;
      iVertex++;
    }
    
    _nRows = nVertices;
    _nCols = nVertices;
  }
  infile.close();
  return allOk;
}

bool AdjSMatGraph::writeEdgeListToFile(std::string filename){
  bool allOk = true;
  std::ofstream outfile;
  
  outfile.open(filename.c_str());
  if(outfile.is_open()){
    for(int iVertex = 0; iVertex < _nRows; iVertex++){
      for(int jVertex = _rowIndex[iVertex]; jVertex < _rowIndex[iVertex+1]; jVertex++){
          outfile << iVertex << " " << _cols[jVertex] << std::endl;
        }
    }
  }else{
    allOk = false;
  }
  return allOk;
}


bool AdjSMatGraph::checkAdjMat(){
  bool allOk = true;

  if(_nRows > 1 && _nCols > 1){
    for(int iRow = 0; iRow < (_nRows/2)+1; iRow++){
      for(int iCol = _rowIndex[iRow]; iCol < _rowIndex[iRow+1]; iCol++){
        bool found = false;
        for(int iCheck = _rowIndex[_cols[iCol]]; iCheck < _rowIndex[_cols[iCol]+1]; iCheck++){
          if(_cols[iCheck] == iRow){
            found = true;
          }
        }
        if(!found){
          std::cout << "Cannot find a matching entry for edge from " << iRow <<  " to " << _cols[iCol] << std::endl;
          allOk = false;
          break;
        }
      }
      if(!allOk){
        break;
      }
    }
  }
  return allOk;
}




void AdjSMatGraph::printAdjMat(){
  int iRow, iCol, iCurrentCol, iData;
  
  iData = 0;
  if(_nRows > 0 && _nCols > 0){
    for(iRow = 0; iRow <_nRows; iRow++){
      std::cout << iRow << " : ";
      iCurrentCol = _cols[iData];
      if(iData != _rowIndex[iRow+1]){
        for(iCol = 0; iCol <_nCols; iCol++){
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
          std::cout << "0 ";
        }
      }
      std::cout << std::endl;
    }
  }
  return;
}




void AdjSMatGraph::printAdjList(){}

int AdjSMatGraph::connectedComps(std::vector<int>& compIndex){
  compIndex.assign(_nRows,-1);
  
  bool unvisited = true;
  int iComponent = 0, iVertex;
  
  iVertex = 0;
  
  while(unvisited){
    dfsStep(compIndex, iVertex, iComponent);
    
    unvisited = false;
    for(iVertex = 1; (iVertex < _nRows) && (unvisited==false); iVertex++){
      
      if(compIndex[iVertex]==-1){
        unvisited = true;
        break;
      }
    }
    iComponent++;
  }
  if(unvisited){
    std::cout << "Something wrong\n";
  }
  return iComponent;
}

void AdjSMatGraph::dfsStep(std::vector<int>& compIndex, int iVertex, int iConnectedComp){
  int iData;
  
  if(compIndex[iVertex]==-1){
    compIndex[iVertex] = iConnectedComp;
    iData = _rowIndex[iVertex];

    while(iData != _rowIndex[iVertex+1]){
      dfsStep(compIndex,_cols[iData],iConnectedComp);
      iData++;
    }
  }
}

Sparse<int> AdjSMatGraph::getAdjSMatInt(){
  std::vector<int>tempdata(_cols.size(),1);
  Sparse<int> adjMat( _nRows, _nRows,tempdata ,_cols, _rowIndex);
  return adjMat;
}

Sparse<double> AdjSMatGraph::getAdjSMatDbl(){
  std::vector<double>tempdata(_cols.size(),1.0);
 
  Sparse<double> adjMat(_nRows, _nRows,  tempdata, _cols, _rowIndex);
  return adjMat;
}


Sparse<int> AdjSMatGraph::getDegreeSMatInt(){
  std::vector<int> degreeArray(_nRows);
  std::vector<int> colsArray(_nRows);
  int iVertex;
  int nDegree;

  nDegree = 0;

  for(iVertex = 0; iVertex < _nRows; iVertex++){
    nDegree = (int)_rowIndex[iVertex+1]- (int)_rowIndex[iVertex];
    degreeArray[iVertex] = nDegree;
  }
  for(iVertex = 0; iVertex < _nRows; iVertex++){
    colsArray[iVertex] = iVertex;
  }
  std::vector<int> indexArray(_cols.size()+1);
  for(int i = 0; i <= (int)_cols.size(); i++){
    indexArray[i] = i;
  }

  Sparse<int> adjMat( _nRows, _nCols, degreeArray,colsArray,indexArray);
  return adjMat;
}

Sparse<double> AdjSMatGraph::getDegreeSMatDbl(){
  std::vector<double> degreeArray(_nRows);
  std::vector<int> colsArray(_nRows);
  int iVertex, nDegree;

 
  for(iVertex = 0; iVertex < _nRows; iVertex++){
    nDegree = _rowIndex[iVertex+1]- _rowIndex[iVertex];
    degreeArray[iVertex] = (double)nDegree;
  }
  for(iVertex = 0; iVertex < _nRows; iVertex++){
    colsArray[iVertex] = iVertex;
  }
  std::vector<int> indexArray(_cols.size()+1);
  for(int i = 0; i <= (int)_cols.size(); i++){
    indexArray[i] = i;
  }

  Sparse<double>adjMat( _nRows, _nCols, degreeArray,colsArray,indexArray);
  
  return adjMat;
}


Sparse<double> AdjSMatGraph::getLaplacian(){
  
  
  std::vector<double>tempData(_nRows+ _cols.size());
  std::vector<int>tempCols(_nRows+ _cols.size());
  std::vector<int>tempRowIndex(_nRows+1);
  int iAdjData = 0, iLaplaceData = 0;
  for(int iRow = 0; iRow < _nRows;iRow++){
    tempRowIndex[iRow] = iLaplaceData;
    iAdjData = _rowIndex[iRow];
    while(iAdjData < _rowIndex[iRow+1] && _cols[iAdjData] < iRow){
      tempData[iLaplaceData] = -1;
      tempCols[iLaplaceData] = _cols[iAdjData];
      iLaplaceData++;
      iAdjData++;
    }
    tempData[iLaplaceData] =  (double)_rowIndex[iRow+1]- _rowIndex[iRow];
    tempCols[iLaplaceData] = iRow;
    iLaplaceData++;
    
    while(iAdjData < _rowIndex[iRow+1]){
      tempData[iLaplaceData] = -1;
      tempCols[iLaplaceData] = _cols[iAdjData];
      iLaplaceData++;
      iAdjData++;
    }
  }
  tempRowIndex[_nRows] = iLaplaceData;
  
  Sparse<double>laplacian(_nRows, _nCols, tempData,tempCols,tempRowIndex);
  return laplacian;
}

