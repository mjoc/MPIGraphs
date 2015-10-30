#ifndef __MOC_DistribGraph_
#define __MOC_DistribGraph_

#include "Graph.hpp"
#include <vector>
#include <sstream>
#include "mpi.h"

class Runlog;
class SparseBinary;

class DistribGraph: public AdjSMatGraph{
protected:

  MPI_Comm _communicator;
  int _mpi_size;
  int _mpi_rank;

  Runlog *_runlog;
  void statusMessage(char* message);
  
  int _nGlobalVertices;
  int _iProcessVertexOffset;
  std::vector<int> _iLocalVertices;
  std::vector<std::vector<int> > _externalEdges;
  std::vector<std::vector<int> > _externalEdgeRanks;
  SparseBinary extendedAdj;
  int _nExternalEdges;
  std::vector<int> _globalVertices;
  
  std::vector<int> _nLocalVerticesByRank;
  std::vector<double> _randomVariables;

  std::vector<int> _graphPartition;
  
  bool readPartitionFromFile(const char* filename);
  void setGraphPartition(std::vector<int> partition);
  
public:
  ~DistribGraph();
  DistribGraph(MPI_Comm communicator, int rank, int size, Runlog* runlog);
  
  
  const std::vector<std::vector<int> > externalEdges() const;
  const std::vector<std::vector<int> > externalEdgeRanks() const;
  
  int nExternalEdges(){return _nExternalEdges;}
  int nGlobalVertices(){return _nGlobalVertices;}
  const std::vector<int>& nLocalVertices() const {return _nLocalVerticesByRank;}
  
  int mpiSize(){return _mpi_size;}
  int mpiRank(){return _mpi_rank;}
  int vertexOffset(){return _iProcessVertexOffset;}
  
  DistribGraph& operator=(const DistribGraph&);
 	 
  bool distributeGraphRankZero(AdjSMatGraph&);
  bool distributeGraphRankZero(AdjSMatGraph& graph, std::vector<int> partition);
  bool distributeGraphRankZero(AdjSMatGraph& graph, const char* partitionfile);
  
  bool receiveGraphRankNonZero();
  void setData();
  
  Sparse<int> getAdjSMatInt();
  Sparse<double> getAdjSMatDbl();
  Sparse<int> getExtendedAdjSMatInt();
  Sparse<double> getExtendedAdjSMatDbl();

  
  Sparse<int> getDegreeSMatInt();
  Sparse<double> getDegreeSMatDbl();
  Sparse<int> getExtendedDegreeSMatInt();
  Sparse<double> getExtendedDegreeSMatDbl();
  
  
  Sparse<double> getLaplacian();
  
  void printExternalEdges();
  
  friend class DistribIsing;
 
};


#endif
