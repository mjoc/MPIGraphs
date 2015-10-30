#ifndef __MOC_Graph_
#define __MOC_Graph_

#include <vector>
#include "Sparse.hpp"
#include "SparseBinary.hpp"

template<typename DataType>
class Sparse;


class Graph{
public:
  virtual ~Graph(){};
  virtual bool readFromFile(const char *filename, int formatcode) = 0;
  virtual void printAdjMat() = 0;
  virtual void printAdjList() = 0;
  virtual int nVertices() = 0;
  virtual int nEdges() = 0;
  virtual int connectedComps(std::vector<int>& compIndex) = 0;
  virtual Sparse<int> getAdjSMatInt() = 0;
  virtual Sparse<int> getDegreeSMatInt() = 0;
  virtual Sparse<double> getAdjSMatDbl() = 0;
  virtual Sparse<double> getDegreeSMatDbl() = 0;
  virtual Sparse<double> getLaplacian() = 0;

  
};


class AdjSMatGraph : public Graph, public SparseBinary{
  void dfsStep(std::vector<int>& compIndex, int iVertex, int iConnectedComp);
public:
  bool readFromFile(const char *filename, int formatcode);
 
  bool readSloppyFromFile(char *filename, int formatCode);
  bool writeEdgeListToFile(std::string filename);
  bool checkAdjMat();
  void printAdjMat();
  void printAdjList();
  int nVertices(){ return _nRows;}
  int nEdges(){ return nData();}
  int connectedComps(std::vector<int>& compIndex);
  Sparse<int> getAdjSMatInt();
  Sparse<int> getDegreeSMatInt();
  Sparse<double> getAdjSMatDbl();
  Sparse<double> getDegreeSMatDbl();
  Sparse<double> getLaplacian();
  
  friend class DistribGraph;
  friend class Ising;
};

#endif
