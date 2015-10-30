#ifndef __MOC_Ising_
#define __MOC_Ising_

#include "Graph.hpp"
#include "RngWrapper.hpp"
#include <fstream>
#include <set>
#include <cmath>
#include <queue>
#include <string>
#include "Runlog.hpp"

#define UNCOLOURED 0
#define GREEN 1
#define RED 2

enum bipartiteState { yes, no, untested};

class Ising{
protected:
  AdjSMatGraph *_graph;
  std::vector<int> _colours;
  std::vector<int> _states;
  std::set<int> _greens;
  std::set<int> _reds;
  bipartiteState _bipartite;
  Runlog  *_runlog;
  
  std::string _name;
public:

  ~Ising(){};
  Ising(AdjSMatGraph* graph, Runlog* runlog, std::string name);
  
  bool initialiseGraphMatrix(int nVertices, int nEdges);
  
  bool isBipartite(int iStartVertex, int iStartColour);
  void printColours();
  bool doIsingSimulation(double JoverKBT, int nBurnIn, int nSamples, int stepSize, RngWrapper& rng);
  const std::set<int>& greens();
  const std::set<int>& reds();
  int nGreen();
  int nRed();
  
};

#endif
