#ifndef __MOC_GRAPH_FACTORY_
#define __MOC_GRAPH_FACTORY_

#include "Graph.hpp"
#include <stdexcept>


class GraphFactory{
  //Graph* createAdjList(){
  //  return new AdjListGraph();
  //}
  //Graph* createAdjMat(){
  //  return new AdjMatGraph();
  //}
  Graph* createAdjSMat(){
    return new AdjSMatGraph();
  }

public:
  Graph* operator()(const int type){
    if(type == 0){
      //return createAdjMat();
    //}else if(type == 1){
    //  return createAdjList();
    }else if(type == 2){
      return createAdjSMat();
    }else{
      throw std::invalid_argument("must be 0 for 'adj mat' or 1 for 'adj list' or 2 for 'adj sparse mat'");
    }
    return NULL;
  }
};
    
#endif
