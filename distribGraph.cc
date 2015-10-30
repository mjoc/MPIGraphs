#include <cstdlib>
#include "DistribGraph.hpp"
#include "Graph.hpp"
#include "Sparse.hpp"
#include "Runlog.hpp"
#include "mpi.h"
#include "RngWrapper.hpp"


DistribGraph::~DistribGraph(){    
}

DistribGraph::DistribGraph( MPI_Comm communicator, int rank, int size, Runlog* runlog){
  
  _communicator = communicator;
  _mpi_rank = rank;
  _mpi_size = size;
  _runlog = runlog;
}

bool DistribGraph::readPartitionFromFile(const char* filename){
  bool allOk = true;
  std::ifstream infile(filename);
  std::string data;
  int nData = 0;
  _runlog->writelog(std::string("Reading from partition file: ")+std::string(filename));
  if(!infile.is_open()){
    allOk = false;
    _runlog->writelog(std::string("Couldn't open partition file!").c_str());
  }else{
    while (! infile.eof() )
    {
     

      getline (infile,data);
      if(data.find_first_not_of(' ') == std::string::npos)
      {
        break;
      }
      _graphPartition.push_back(atoi(data.c_str()));
      nData++;
    }
    _runlog->writelog(std::string("Read partition file"));
    infile.close();
  }
  
  return allOk;
  
}

void DistribGraph::setGraphPartition(std::vector<int> partition){
  _graphPartition = partition;
}

const std::vector<std::vector<int> > DistribGraph::externalEdges() const {
  return _externalEdges;
}

const std::vector<std::vector<int> > DistribGraph::externalEdgeRanks() const {
  return _externalEdgeRanks;
}

DistribGraph& DistribGraph::operator=(const DistribGraph&){
  return *this;
}

void DistribGraph::statusMessage(char* message){
  if(_runlog != NULL){
    _runlog->writelog(message);
  }else{
    std::cout << message << std::endl;
  }
} 

bool DistribGraph::distributeGraphRankZero(AdjSMatGraph& graph, std::vector<int> partition){
  bool allOk;
  if((int)partition.size() == graph.nVertices()){
    setGraphPartition(partition);
    allOk = distributeGraphRankZero(graph);
  }else{
    allOk = false;
  }
  return allOk;
}

bool DistribGraph::distributeGraphRankZero(AdjSMatGraph& graph, const char* partitionFile){
  bool allOk;

  allOk = readPartitionFromFile(partitionFile);
  if(allOk){
    allOk = distributeGraphRankZero(graph);
  }
  return allOk;
}


bool DistribGraph::distributeGraphRankZero(AdjSMatGraph& graph){
  bool allOk = true;
  char message[200];
 
  ////////////////////////////////////////////
  bool partitionSupplied = false;
  if(_graphPartition.size()>0){
    if((int)_graphPartition.size() == graph.nVertices()){
      std::string message = ("Using supplied graph partition");
      _runlog->writelog(message);
      partitionSupplied = true;
    }else{
      std::cerr << "Supplied graph partition is not correct (wrong size)\n";
      std::cerr << "#Vertices: " << graph.nVertices() << " Partition Size: " << _graphPartition.size() << "\n";
      allOk = false;
    }
  }
  
  RngWrapper randomGenerator;
  std::vector<int> nLocalDataLengths(_mpi_size);
  std::vector<int> newVertexOrder;
  
  _nLocalVerticesByRank.resize(_mpi_size);
  
  if(partitionSupplied){
    sprintf (message, "Graph partition supplied");
    _runlog->writelog(message);
    newVertexOrder.resize(graph.nVertices());
    std::fill(newVertexOrder.begin(),newVertexOrder.end(),-1);
    
    std::vector<std::vector<int> >vertexByRank(_mpi_size);
    for(int iVertex = 0; iVertex < graph.nVertices();iVertex++){
      if(_graphPartition[iVertex] >= _mpi_size){
        allOk = false;
        std::cout << "Partion Size greater than MPI Processes! \n";
        break;
      }
      vertexByRank[_graphPartition[iVertex]].push_back(iVertex);
    }
    if(allOk){
      int index=0;
      for(int iRank = 0; iRank < _mpi_size; iRank++){
        for(int iVertexInRank = 0; iVertexInRank < (int)vertexByRank[iRank].size(); iVertexInRank++){
          
          newVertexOrder[index] = vertexByRank[iRank][iVertexInRank];
          index++;
        }
        _nLocalVerticesByRank[iRank] = (int)vertexByRank[iRank].size();
      }
    }
  
    for(int i = 0 ; i < (int)newVertexOrder.size(); i++){
      if(newVertexOrder[i] == -1){
         std::cerr << "Problem with graph partition, not all vertices assigned, giving up!\n";
        allOk = false;
      }
    }
  }else{
      sprintf (message, "No graph partition supplied, partitioning by vertex order");
    _runlog->writelog(message);

    newVertexOrder.resize(graph.nVertices());
    for(int i = 0; i < (int)newVertexOrder.size();i++){
      newVertexOrder[i] = i;
    }
    
    //Calculate the number of vertices per process
    int baseSize = (int)graph.nVertices()/_mpi_size;
    int residual = (int)graph.nVertices() - _mpi_size * baseSize;
    {
    int i;
    for(i = 0; i < _mpi_size; i++){
      _nLocalVerticesByRank[i] = baseSize;
    }
    
    i = 0;
    while(residual > 0 && i  < _mpi_size){
      _nLocalVerticesByRank[i]++;
      i++;
      residual--;
    }
    }
  }
  
  std::vector<int> newVertexSwitchOrder(newVertexOrder.size());
  for(int i = 0; i < (int)newVertexOrder.size(); i++){
    newVertexSwitchOrder[newVertexOrder[i]]= i;
  }
  
  int iOffset = 0;
  for(int i = 0; i < _mpi_size; i++){
    for(int j = iOffset; j < iOffset+ _nLocalVerticesByRank[i]; j++){
      int jPermuted = newVertexOrder[j];

      nLocalDataLengths[i] += (int)(graph._rowIndex[jPermuted+1] - graph._rowIndex[jPermuted]);
    }
    iOffset += _nLocalVerticesByRank[i];
  }
  
  _nGlobalVertices = graph.nVertices();

  sprintf (message, "Distribution of Vertices by Rank");
  _runlog->writelog(message);
  for(size_t i = 0; i < _nLocalVerticesByRank.size();i++){
    sprintf (message, "%lu %d",i,_nLocalVerticesByRank[i]);
  _runlog->writelog(message);
  }

  
  int nLocalData;
  
  MPI_Bcast(&_nGlobalVertices, 1, MPI_INT, 0, _communicator);
  MPI_Bcast(&_nLocalVerticesByRank[0], _mpi_size, MPI_INT, 0,_communicator);
  MPI_Scatter(&nLocalDataLengths[0], 1, MPI_INT, &nLocalData,1,MPI_INT,0,_communicator);
 
 
  std::vector<int> newCols;
  std::vector<int> newRowIndex;
  std::vector<int> originalVertex;
  int iVertexOffset = _nLocalVerticesByRank[0];;
  int jPermuted;

  for(int iRank = 1; iRank < _mpi_size; iRank++){
    newCols.resize(nLocalDataLengths[iRank]);
    newRowIndex.resize(_nLocalVerticesByRank[iRank]+1);
    
    originalVertex.resize(_nLocalVerticesByRank[iRank]);

    int iIndex = 0;
    for(int iVertex = iVertexOffset; iVertex < iVertexOffset + _nLocalVerticesByRank[iRank]; iVertex++){
      jPermuted = newVertexOrder[iVertex];
      
      newRowIndex[iVertex-iVertexOffset] = iIndex;
      originalVertex[iVertex - iVertexOffset] = jPermuted;
      std::set<int> orderedCols;
      for(int iRowData = graph._rowIndex[jPermuted]; iRowData< graph._rowIndex[jPermuted+1]; iRowData++){
        orderedCols.insert(newVertexSwitchOrder[graph._cols[iRowData]]);
      }
      //sprintf (message, "Adding %ld to adjancy matrix: ",orderedCols.size());
      //_runlog->writelog(message);
      for(std::set<int>::const_iterator it = orderedCols.begin(); it != orderedCols.end(); ++it){
        newCols[iIndex] = *it;
        iIndex++;
      }
    }
    newRowIndex[_nLocalVerticesByRank[iRank]] = iIndex;
    iVertexOffset += _nLocalVerticesByRank[iRank];

    MPI_Send(&newCols[0],nLocalDataLengths[iRank],MPI_INT, iRank, 1, _communicator);
    MPI_Send(&newRowIndex[0],(int)_nLocalVerticesByRank[iRank]+1,MPI_INT, iRank, 2, _communicator);
    MPI_Send(&originalVertex[0],(int)_nLocalVerticesByRank[iRank],MPI_INT, iRank, 3, _communicator);
  }
  
  extendedAdj._cols.resize(nLocalData);
  extendedAdj._rowIndex.resize(_nLocalVerticesByRank[0]+1);
  _globalVertices.resize(_nLocalVerticesByRank[0]);
  int iIndex = 0;
  for(int iVertex = 0; iVertex < _nLocalVerticesByRank[0]; iVertex++){
    jPermuted = newVertexOrder[iVertex];
    extendedAdj._rowIndex[iVertex] = iIndex;
    _globalVertices[iVertex] = jPermuted;
    std::set<int> orderedCols;
    for(int iRowData = graph._rowIndex[jPermuted]; iRowData< graph._rowIndex[jPermuted+1]; iRowData++){
      orderedCols.insert(newVertexSwitchOrder[graph._cols[iRowData]]);
    }
    for(std::set<int>::const_iterator it = orderedCols.begin(); it != orderedCols.end(); ++it){
      extendedAdj._cols[iIndex] = *it;
      iIndex++;
    }
  }
  extendedAdj._rowIndex[_nLocalVerticesByRank[0]] = iIndex;
  
  setData();
  
  return true;
}

bool DistribGraph::receiveGraphRankNonZero(){
  int nLocalVertices = 0, dummy = 0;
  int dummy2 = 0;
  int nLocalData;
  MPI_Status status;
  char message[200];
  
  MPI_Scatter(&dummy, dummy2, MPI_INT, &nLocalData,1,MPI_INT,0,_communicator);
  
  
  _nLocalVerticesByRank.resize(_mpi_size);

  MPI_Bcast(&_nGlobalVertices, 1, MPI_INT, 0,_communicator);
  sprintf(message, "Number of Global Vertices is: %d", _nGlobalVertices);
  _runlog->writelog(message);

  MPI_Bcast(&_nLocalVerticesByRank[0], _mpi_size, MPI_INT, 0,_communicator);
  sprintf(message, "Number of Local Vertices is: %d", _nLocalVerticesByRank[_mpi_rank]);
  _runlog->writelog(message);
  
  //MPI_Scatter(&dummy, dummy, MPI_INT, &nLocalData,1,MPI_INT,0,_communicator);

  sprintf(message, "Number of local data is: %d", nLocalData);
  _runlog->writelog(message);

  
  nLocalVertices = _nLocalVerticesByRank[_mpi_rank];

  extendedAdj._cols.resize(nLocalData);
  extendedAdj._rowIndex.resize(nLocalVertices+1);
  _globalVertices.resize(nLocalVertices);
  
  ///
  int nReceive = nLocalData;

  MPI_Recv(&extendedAdj._cols[0], nReceive, MPI_INT, 0, 1, _communicator, &status);

  MPI_Recv(&extendedAdj._rowIndex[0], nReceive, MPI_INT, 0, 2, _communicator, &status);

  MPI_Recv(&_globalVertices[0], nReceive, MPI_INT, 0, 3, _communicator, &status);
  
  setData();

  return true;
}

void DistribGraph::setData(){
  _iProcessVertexOffset = 0;
  {
    int i = 0;
    while(i < _mpi_rank){
      _iProcessVertexOffset += _nLocalVerticesByRank[i];
      i++;
    }
  }

  int nEdgesInAndEx = (int)extendedAdj._cols.size();
  int nEdges = 0;
  for(int i = 0; i < nEdgesInAndEx; i++){
    if(extendedAdj._cols[i] >= _iProcessVertexOffset &&  extendedAdj._cols[i] < _iProcessVertexOffset + _nLocalVerticesByRank[_mpi_rank]){
      nEdges++;
    }
  }

  int nLocalVertices = _nLocalVerticesByRank[_mpi_rank];

  _externalEdges.resize(nLocalVertices);
  _externalEdgeRanks.resize(nLocalVertices);
  _nRows = nLocalVertices;
  _nCols = nLocalVertices;
  _cols.resize(nEdges);
  _rowIndex.resize(_nRows+1);

  extendedAdj._nRows = nLocalVertices;
  extendedAdj._nCols = _nGlobalVertices ;
 
   std::string messageString = std::string("There are: ") + _runlog->to_string(_nRows) + std::string(" local vertices of  ") + _runlog->to_string(extendedAdj._nCols) + " total";
  _runlog->writelog(messageString);
  
  
  //Figure out where these external vertices are
  //std::vector<int> nLocalVertices = graph->getNLocalVertices();
  std::vector<int> nCummulativeLocalVertices(_mpi_size);
  nCummulativeLocalVertices[0] = _nLocalVerticesByRank[0];
  for(int i = 1; i < _mpi_size; i++){
    nCummulativeLocalVertices[i] = nCummulativeLocalVertices[i-1]+ _nLocalVerticesByRank[i];
  }
  
  _nExternalEdges = 0;
  _rowIndex[0] = 0;
  extendedAdj._rowIndex[0] = 0;

  int iIn = 0, iNextRow = 1;
  int iRank;
  for(int iEdge = 0; iEdge < nEdgesInAndEx; iEdge++){
    while(iEdge >= extendedAdj._rowIndex[iNextRow]){
      _rowIndex[iNextRow] = iIn;
      iNextRow++;
    }

    if(extendedAdj._cols[iEdge] >= _iProcessVertexOffset &&  extendedAdj._cols[iEdge] < _iProcessVertexOffset+ nLocalVertices){
      _cols[iIn] = extendedAdj._cols[iEdge]- _iProcessVertexOffset;
       ++iIn;
    }else{
      _externalEdges[iNextRow-1].push_back(extendedAdj._cols[iEdge]);
      iRank = 0;
      while(extendedAdj._cols[iEdge] >= nCummulativeLocalVertices[iRank]){
        iRank++;
      }
      _externalEdgeRanks[iNextRow-1].push_back(iRank);
      _nExternalEdges++;
    }
  }
  messageString = std::string("There are: ") + _runlog->to_string(_nExternalEdges) + " external edges";
  _runlog->writelog(messageString);

  std::vector<int> externalEdgesAcrossRanks(_mpi_size);
  MPI_Gather(&_nExternalEdges,1,MPI_INT,&externalEdgesAcrossRanks[0],1,MPI_INT,0,_communicator);

  _runlog->writelog(std::string("External Edges per rank:").c_str());

  for(int i = 0;i < _mpi_size; i++){
    messageString = _runlog->to_string(i) + " " + _runlog->to_string(externalEdgesAcrossRanks[i]);
    _runlog->writelog(messageString);
  }
  //int MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
  //             void *recvbuf, int recvcount, MPI_Datatype recvtype,
  //             int root, MPI_Comm comm)

  //iIn++;
  while(iNextRow <= _nRows){
    _rowIndex[iNextRow] =iIn;
    iNextRow++;
  }
   
  return;
}

void DistribGraph::printExternalEdges(){
  char message[200];
  snprintf (message,200, "There are %d external edges",_nExternalEdges);
  _runlog->writelog(message); 

  for(int i = 0; i < (int)_externalEdges.size(); i++){
    snprintf (message,200, "Vertex %d: ",i);
    for(std::vector<int>::iterator it = _externalEdges[i].begin(); it != _externalEdges[i].end();++it){
      snprintf (message,200, "%s %d",message,(int)*it);
    }
    _runlog->writelog(message);
    
  }

  sprintf (message, "Ranks: ");
  _runlog->writelog(message);
  
  for(int i = 0; i < (int)_externalEdgeRanks.size(); i++){
    sprintf (message, "Vertex %d: ",i);
    for(std::vector<int>::iterator it = _externalEdgeRanks[i].begin(); it != _externalEdgeRanks[i].end();++it){
      sprintf (message, "%s %d",message,(int)*it);
    }
    _runlog->writelog(message);
    
  }
}

Sparse<int> DistribGraph::getAdjSMatInt(){
  std::vector<int>tempdata(_cols.size(),1);
  Sparse<int> adjMat( _nRows, _nCols,  tempdata, _cols, _rowIndex);
  return adjMat;
}

Sparse<int> DistribGraph::getExtendedAdjSMatInt(){
  std::vector<int>tempdata(extendedAdj._nCols,1);
  Sparse<int> adjMat( _nRows, extendedAdj._nCols, tempdata, extendedAdj._cols, extendedAdj._rowIndex);
  return adjMat;
}

Sparse<double> DistribGraph::getAdjSMatDbl(){
  std::vector<double>tempdata(_cols.size(),1.0);

  Sparse<double> adjMat(_nRows, _nCols, tempdata, _cols, _rowIndex);
  return adjMat;
}

Sparse<double> DistribGraph::getExtendedAdjSMatDbl(){
  std::vector<double>tempdata(extendedAdj._nCols,1.0);
   
  Sparse<double> adjMat(_nRows, extendedAdj._nCols, tempdata, extendedAdj._cols, extendedAdj._rowIndex);
  return adjMat;
}

Sparse<int> DistribGraph::getDegreeSMatInt(){
  std::vector<int> degreeArray(_nRows);
  std::vector<int> colsArray(_nRows);
  int iVertex = 0,  nDegree = 0;
  
  for(iVertex = 0; iVertex < _nRows; iVertex++){
    nDegree =  extendedAdj._rowIndex[iVertex+1] -  extendedAdj._rowIndex[iVertex];
    degreeArray[iVertex] = (int)nDegree;
    colsArray[iVertex] = iVertex;
  }
 
  std::vector<int> indexArray(_nRows+1);
  iVertex = 0;
  for(std::vector<int>::iterator it = indexArray.begin(); it != indexArray.end(); ++it){
    *it = iVertex;
    iVertex++;
  }
  
  Sparse<int>adjMat( _nRows, _nCols, degreeArray, colsArray, indexArray);
  return adjMat;
}


Sparse<int> DistribGraph::getExtendedDegreeSMatInt(){
  std::vector<int> degreeArray(_nRows);
  std::vector<int> colsArray(_nRows);
  int iVertex = 0;
  
  for(iVertex = 0; iVertex < _nRows; iVertex++){
    degreeArray[iVertex] = (int)(extendedAdj._rowIndex[iVertex+1] -  extendedAdj._rowIndex[iVertex]);
    colsArray[iVertex] = iVertex+ _iProcessVertexOffset;
  }

  std::vector<int> indexArray(_nRows+1);
  iVertex = 0;
  for(std::vector<int>::iterator it = indexArray.begin(); it != indexArray.end(); ++it){
    *it = iVertex;
    iVertex++;
  }
  
  Sparse<int>adjMat( _nRows, extendedAdj._nCols, degreeArray,colsArray,indexArray);
  return adjMat;
}

Sparse<double> DistribGraph::getDegreeSMatDbl(){
  std::vector<double> degreeArray(_nRows);
  std::vector<int> colsArray(_nRows);
  int iVertex = 0;
  
  for(iVertex = 0; iVertex < _nRows; iVertex++){
    degreeArray[iVertex] = (double)extendedAdj._rowIndex[iVertex+1] -  extendedAdj._rowIndex[iVertex];
    colsArray[iVertex] = iVertex;
  }

  std::vector<int> indexArray(_nRows+1);
  iVertex = 0;
  for(std::vector<int>::iterator it = indexArray.begin(); it !=  indexArray.end(); ++it){
    *it = iVertex;
    iVertex++;
  }
  
  Sparse<double>adjMat( _nRows, _nRows, degreeArray,colsArray,indexArray);
  return adjMat;
}

Sparse<double> DistribGraph::getExtendedDegreeSMatDbl(){
  std::vector<double> degreeArray(_nRows);
  std::vector<int> colsArray(_nRows);
  int iVertex = 0;
  
  for(iVertex = 0; iVertex < _nRows; iVertex++){
    degreeArray[iVertex] = (double)extendedAdj._rowIndex[iVertex+1] -  extendedAdj._rowIndex[iVertex];
    colsArray[iVertex] = iVertex+ _iProcessVertexOffset;
  }

  std::vector<int> indexArray(_nRows+1);
  iVertex = 0;
  for(std::vector<int>::iterator it = indexArray.begin(); it != indexArray.end(); ++it){
    *it = iVertex;
    iVertex++;
  }
  
  Sparse<double>adjMat( _nRows,extendedAdj._nCols, degreeArray,colsArray,indexArray);
  return adjMat;
}



Sparse<double> DistribGraph::getLaplacian(){
  std::vector<double>tempData(extendedAdj._nRows+ extendedAdj._cols.size());
  std::vector<int>tempCols(extendedAdj._nRows+ extendedAdj._cols.size());
  std::vector<int>tempRowIndex(extendedAdj._nRows+1);
  
  int iAdjData = 0, iLaplaceData = 0;
  for(int iRow = 0; iRow < extendedAdj._nRows;iRow++){
    tempRowIndex[iRow] = iLaplaceData;
    iAdjData = extendedAdj._rowIndex[iRow];
    while(iAdjData < extendedAdj._rowIndex[iRow+1] && extendedAdj._cols[iAdjData] < iRow){
      tempData[iLaplaceData] = -1;
      tempCols[iLaplaceData] = extendedAdj._cols[iAdjData];
      iLaplaceData++;
      iAdjData++;
    }
    tempData[iLaplaceData] =  (double)extendedAdj._rowIndex[iRow+1]- extendedAdj._rowIndex[iRow];
    tempCols[iLaplaceData] = iRow;
    iLaplaceData++;
    
    while(iAdjData < extendedAdj._rowIndex[iRow+1]){
      tempData[iLaplaceData] = -1;
      tempCols[iLaplaceData] = extendedAdj._cols[iAdjData];
      iLaplaceData++;
      iAdjData++;
    }
  }
  tempRowIndex[_nRows] = iLaplaceData;
  
  Sparse<double>laplacian(extendedAdj._nRows, extendedAdj._nCols, tempData,tempCols,tempRowIndex);
  return laplacian;
}

