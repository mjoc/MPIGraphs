#ifndef __MOC_Distrib_Ising_
#define __MOC_Distrib_Ising_

#include "DistribGraph.hpp"
#include "Ising.hpp"
#include <vector>
#include "mpi.h"

template<typename DataType>
class My_sparse;

class Runlog;

class DistribGraph;

#define UNCOLOURED 0
#define GREEN 1
#define RED 2


class DistribIsing: public Ising{
  MPI_Comm _communicator;
  Runlog *_runlog;
  DistribGraph *_distribGraph;
  int _mpi_rank;
  int _mpi_size;
  bool _mpi_ok;
  std::string _directory;

  std::vector<int> _sendDataInVertex;
  
  std::vector<int> _sendDataInVertexGlobal;
  std::vector<int> _sendDataExVertexGlobal;
  std::vector<int> _sendDataExVertexRank;
  std::vector<int> _sendData;
  std::vector<int> _recvData;
  std::vector<int> _sendDataNByRank;
  std::vector<int> _recvDataNByRank;
  
  std::vector<int> _sendDataRankDisplacements;
  std::vector<int> _recvDataRankDisplacements;
  
  std::vector<std::vector<int> > _recvDataInVertices;
  std::vector<std::vector<int> > _recvDataIndexByInVertex;
  
  
  std::vector<int> _sendDataRed;
  std::vector<int> _recvDataRed;
  std::vector<int> _sendDataGreen;
  std::vector<int> _recvDataGreen;
  
  std::vector<int> _sendDataGreenInVertex;
  std::vector<int> _sendDataGreenInVertexGlobal;
  std::vector<int> _sendDataRedInVertex;
  std::vector<int> _sendDataRedInVertexGlobal;
  
  std::vector<int> _sendDataGreenNByRank;
  std::vector<int> _recvDataGreenNByRank;
  std::vector<int> _sendDataRedNByRank;
  std::vector<int> _recvDataRedNByRank;

  std::vector<int> _sendDataGreenRankDisplacements;
  std::vector<int> _recvDataGreenRankDisplacements;
  std::vector<int> _sendDataRedRankDisplacements;
  std::vector<int> _recvDataRedRankDisplacements;

  std::vector<int> _recvDataGreenToFullIndex;
  std::vector<int> _recvDataRedToFullIndex;
  double _timer;
  double _timerPackUnpack;
  double _timerAll2All;
  double _dataVolumeIn;
  double _dataVolumeOut;
public:

  DistribIsing( int rank, int size, MPI_Comm comm, DistribGraph* graph, Runlog *runlog, std::string modelName ,std::string fileDirectory):
  Ising((AdjSMatGraph*) graph, runlog, modelName){
    //char message[200];
    
    _distribGraph = graph;
    _mpi_rank = rank;
    _mpi_size = size;
    _communicator = comm;
    _runlog = runlog;
    _directory = fileDirectory;
    
    const std::vector<std::vector<int> > externalEdges = graph->externalEdges();
    const std::vector<std::vector<int> > externalEdgeRanks = graph->externalEdgeRanks();
    
    int iIndex =0, nRecvData = 0, nSendData = 0;
    std::set<int> externalEdgeUniqueInVertexByRank;
    std::map<int,int> externalEdgeUniqueExVertexByRank;

    for(int iRank = 0; iRank < size; iRank++){
      externalEdgeUniqueInVertexByRank.clear();
      externalEdgeUniqueExVertexByRank.clear();
      for(int iInternal = 0; iInternal < (int)externalEdgeRanks.size(); iInternal++){
        for(int iExternal = 0; iExternal < (int)externalEdgeRanks[iInternal].size();iExternal++){

          if(externalEdgeRanks[iInternal][iExternal]== iRank){
            externalEdgeUniqueExVertexByRank.insert(std::make_pair(externalEdges[iInternal][iExternal],iRank));
            externalEdgeUniqueInVertexByRank.insert(iInternal);
            iIndex++;
          }
        }
      }
      nRecvData += (int)externalEdgeUniqueExVertexByRank.size();
      nSendData += (int)externalEdgeUniqueInVertexByRank.size();
      
      _recvDataNByRank.push_back((int)externalEdgeUniqueExVertexByRank.size());
      for(std::map<int,int>::const_iterator it = externalEdgeUniqueExVertexByRank.begin(); it != externalEdgeUniqueExVertexByRank.end(); it++){
        _sendDataExVertexGlobal.push_back(it->first);
        _sendDataExVertexRank.push_back(it->second);
      }
      
      _sendDataNByRank.push_back((int)externalEdgeUniqueInVertexByRank.size());
      
      for(std::set<int>::const_iterator it = externalEdgeUniqueInVertexByRank.begin(); it != externalEdgeUniqueInVertexByRank.end(); it++){
        _sendDataInVertexGlobal.push_back(graph->_globalVertices[*it]);
        _sendDataInVertex.push_back(*it);
      }
    }
   
    _sendData.resize(nSendData);
    _recvData.resize(nRecvData);
    
    _sendDataRankDisplacements.resize(size);
    for(int iRank = 1; iRank < size; iRank++){
      _sendDataRankDisplacements[iRank] = _sendDataRankDisplacements[iRank-1] + _sendDataNByRank[iRank-1];
    }
    _recvDataRankDisplacements.resize(size);
    for(int iRank = 1; iRank < size; iRank++){
      _recvDataRankDisplacements[iRank] = _recvDataRankDisplacements[iRank-1] + _recvDataNByRank[iRank-1];
    }
   
    _recvDataInVertices.resize(nRecvData);
    
    for(int iExVertex = 0; iExVertex < (int)_sendDataExVertexGlobal.size(); iExVertex++){
      for(int iInternal = 0; iInternal < (int)externalEdges.size(); iInternal++){
        for(int iExternal = 0; iExternal < (int)externalEdges[iInternal].size();iExternal++){
          if(externalEdges[iInternal][iExternal] == _sendDataExVertexGlobal[iExVertex]){
            _recvDataInVertices[iExVertex].push_back(iInternal);
          }
        }
      }
    }
    
    _recvDataIndexByInVertex.resize(graph->nVertices());
    {
      int i = 0;
      for(std::vector<std::vector<int> >::const_iterator it = _recvDataInVertices.begin(); it != _recvDataInVertices.end(); ++it){
        for(std::vector<int>::const_iterator it2 = (*it).begin(); it2!= (*it).end(); ++it2){
          _recvDataIndexByInVertex[*it2].push_back(i);
        }
        i++;
      }
    }
      
    
    MPI_Alltoallv(&_sendDataInVertexGlobal[0],&_sendDataNByRank[0],&_sendDataRankDisplacements[0], MPI_INT,&_recvData[0],&_recvDataNByRank[0],&_recvDataRankDisplacements[0], MPI_INT, _communicator);
    
  
  };

  ~DistribIsing(){}
  void statusMessage(char* message);
  bool isBipartite(int iStartVertex);
  bool doIsingSimulation(double JoverKBT, int nBurnIn, int nSamples, int stepSize, RngWrapper& rng);
  void sendAndReceive(std::vector<int> sendVector);
  void sendAndReceive(std::vector<int> sendVector, int colour);
  void printReceiveDataInternalVertices();
};

#endif
