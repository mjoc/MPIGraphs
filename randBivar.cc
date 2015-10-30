#include <fstream>
#include <iostream>
#include <cmath>
#include <set>
#include <algorithm>
#include "RngWrapper.hpp"

#define nGroups 16
#define groupIntConnLow 4
#define groupIntConnHigh 10

#define groupExtConnLow 6
#define groupExtConnHigh 10

//#define groupExtConnLowProb 0.01
//#define groupExtConnHighProb 0.01

#define groupSizeLow 400
#define groupSizeHigh 500


int main(int argc, char **argv){
  RngWrapper rands;
  rands.setSeed(12039);
  std::vector<double> runif;
  int nVertices = 0;
  std::vector<int>nGroupSizes(nGroups);
  std::vector<int>iGroupStarts(nGroups+1);
  if(groupSizeLow == groupSizeHigh){
    for(int i = 0; i < nGroupSizes.size();i++){
      nGroupSizes[i] = groupSizeLow;
      nVertices+= nGroupSizes[i];
      if(i == 0){
        iGroupStarts[i] = 0;
      }else{
        iGroupStarts[i] = iGroupStarts[i-1] + nGroupSizes[i-1];
      }
      
    }
  }else{
    runif = rands.getUniform01Randoms(nGroups);
    for(int i = 0; i < nGroupSizes.size();i++){
      nGroupSizes[i] = groupSizeLow + runif[i]*(groupSizeHigh - groupSizeLow + 1) ;
      nVertices+= nGroupSizes[i];
      if(i == 0){
        iGroupStarts[i] = 0;
      }else{
        iGroupStarts[i] = iGroupStarts[i-1] + nGroupSizes[i-1];
      }
    }
    iGroupStarts[nGroups] = nVertices;
  }
  /*
  std::vector<double>groupExtEdgeProb(nGroups);
  if(groupExtConnLowProb == groupExtConnHighProb){
    for(int i = 0; i < nGroupSizes.size();i++){
      groupExtEdgeProb[i] = groupExtConnLowProb;
    }
  }else{
    runif = rands.getUniform01Randoms(nGroups);
    for(int i = 0; i < nGroupSizes.size();i++){
      groupExtEdgeProb[i] = groupExtConnLowProb + runif[i]*(groupExtConnHighProb - groupExtConnLowProb) ;
    }
  }
  */
  int nGreen = 0, nRed = 0;
  runif = rands.getUniform01Randoms(nVertices);
  std::vector<int> redGreen(nVertices);
  std::vector<std::vector<int> > greenByGroup(nGroups);
  std::vector<std::vector<int> > redByGroup(nGroups);
  
  for(int iGroup = 0; iGroup < nGroups;iGroup++){
    for(int iVertex = iGroupStarts[iGroup]; iVertex < iGroupStarts[iGroup+1]; iVertex++){
      if(runif[iVertex]>0.5){
        redGreen[iVertex] = 0;
        greenByGroup[iGroup].push_back(iVertex);
        nGreen++;
      }else{
        redGreen[iVertex] = 1;
        redByGroup[iGroup].push_back(iVertex);
        nRed++;
      }
    }
  }
  std::cout << "Green: " << nGreen << " Red: " << nRed << " total: " << nVertices  <<std::endl;
  
  std::set< std::pair<int,int> > edges;
  int nInt = 0, nExt = 0;
  
  for(int iGroup1 = 0; iGroup1 < nGroups; iGroup1++){
    
    int nConnects = groupIntConnLow + (int)floor(rands.getUniform01Random()*(groupIntConnHigh-groupIntConnLow+1));
    //std::cout << "Group: " << iGroup1 << " internal connections: " << nConnects << "\n";
    for(int iVertex1 = iGroupStarts[iGroup1]; iVertex1 < iGroupStarts[iGroup1+1];iVertex1++){
      
      if(redGreen[iVertex1]==0){
        //nConnects = std::min(2,3);
        nConnects = std::min((int)redByGroup[iGroup1].size(),nConnects );
      }else{
        nConnects = std::min((int)greenByGroup[iGroup1].size(),nConnects );
      }
      
      std::vector<int>iConnected =  rands.getRandomPermutation(nConnects);
      int iVertex2;
      for(int i = 0; i < nConnects; i++){
        if(redGreen[iVertex1]==0){
          iVertex2 = redByGroup[iGroup1][iConnected[i]];
        }else{
          iVertex2 = greenByGroup[iGroup1][iConnected[i]];
        }
        if(iVertex1 < iVertex2){
          if(edges.insert(std::make_pair(iVertex1,iVertex2)).second){
            nInt++;
            //std::cout << "  Group: " << iGroup1 <<" internal : " << iVertex1 << " - " << iVertex2 << "\n:";
          }
        }else{
          if(edges.insert(std::make_pair(iVertex2,iVertex1)).second){
            nInt++;
            //std::cout << "  Group: " << iGroup1 <<" internal : " << iVertex2 << " - " << iVertex1 << "\n:";
          }
        }
      }
    }
    
    nConnects = groupIntConnLow + (int)floor(rands.getUniform01Random()*(groupExtConnHigh-groupExtConnLow+1));
    nConnects = std::min(nConnects, nGroupSizes[iGroup1]);
    std::vector<int> extGroupInt = rands.getRandomPermutation(nConnects);
    
    std::vector<double> randoms =  rands.getUniform01Randoms(nConnects);
    std::vector<int> extGroup(randoms.size());
    for(int i = 0; i < randoms.size();i++){
      extGroup[i] = (int)floor(randoms[i]* (nGroups-1));
      if(extGroup[i] >= iGroup1){
        extGroup[i]++;
      }
    }
    
    //std::cout << "Group: " << iGroup1 << " external connections: " << nConnects << "\n";
    
    for(int iExEdges = 0; iExEdges < extGroupInt.size();iExEdges++){
      int iVertex1 = extGroupInt[iExEdges] + iGroupStarts[iGroup1];
      int iGroup2 = extGroup[iExEdges];
      int iVertex2 = -1;
      if(redGreen[iVertex1] == 0){
        iVertex2 = redByGroup[iGroup2][(int)floor(redByGroup[iGroup2].size()*rands.getUniform01Random())];
      }else{
        iVertex2 = greenByGroup[iGroup2][(int)floor(greenByGroup[iGroup2].size()*rands.getUniform01Random())];
      }
      if(iVertex1 < iVertex2){
        if(edges.insert(std::make_pair(iVertex1,iVertex2)).second){
          nExt++;
          //std::cout << "  Group: " << iGroup1 <<" to external : " << iGroup2 << " Vertex: " << iVertex1 << " - " << iVertex2 << "\n:";
        }
      }else{
        if(edges.insert(std::make_pair(iVertex2,iVertex1)).second){
          nExt++;
          //std::cout << "  Group: " << iGroup1 <<" to external : " << iGroup2 << " Vertex: " << iVertex2 << " - " << iVertex1 << "\n:";
        }
      }
    }
  }
  
  std::cout << "Total Edges: " << edges.size() << " on " << nVertices << " vertices\n";
  std::cout << "Internal: " << nInt << " external: " << nExt << "\n";

  char nVerticesChar[20];
  sprintf(nVerticesChar, "%d", nVertices);
  
  char nEdgesChar[20];
  sprintf(nEdgesChar, "%d", nInt + nExt);
  
  char nGroupsChar[20];
  sprintf(nGroupsChar, "%d", nGroups);
  
  
  std::string outDir("/home/users/mschpc/2013/oconnm28/Courses/GraphPartitioning/Code/Graphs/");
  std::string fileName(outDir + std::string("bivariateV") + nVerticesChar + "E" + nEdgesChar + "G" + nGroupsChar + "_el.txt");
  std::ofstream elFile;
  elFile.open(fileName.c_str());

  if(!elFile.is_open()){
    std::cerr<< "problem opening el file, giving up!\n";
    exit(1);
  }

  std::vector<std::vector<int> > adjMat(nVertices);
  for(int i = 0; i < nVertices; i++){
    adjMat[i].resize(nVertices);
  }
  
  for(std::set<std::pair<int, int> >::iterator it = edges.begin(); it !=  edges.end(); ++it){
    elFile << (*it).first << " " << (*it).second << "\n";
    if((*it).second <= (*it).first){
      std::cerr << "Something wrong\n";
      exit(1);
    }else{
      adjMat[(*it).first][(*it).second] = 1;
      adjMat[(*it).second][(*it).first] = 1;
    }
  }
  elFile.close();
  
  
  
  fileName = outDir + std::string("bivariateV") + nVerticesChar + "E" + nEdgesChar + "G" + nGroupsChar + ".txt";
  
  std::ofstream adjFile;
  adjFile.open(fileName.c_str());
  if(!adjFile.is_open()){
    std::cerr<< "problem opening adj file, giving up!\n";
    exit(1);
  }

  for(int i = 0; i < nVertices; i++){
    adjFile << i ;
    for(int j = 0; j < nVertices; j++){
      if(adjMat[i][j] == 1){
        adjFile << "; " << j;
      }
    }
    adjFile << std::endl;
  }
  adjFile.close();
  
  fileName = outDir + std::string("bivariateV") + nVerticesChar + "E" + nEdgesChar + "G" + nGroupsChar + "_gps.txt";
  
  std::ofstream groupsFile;
  groupsFile.open(fileName.c_str());
  if(!groupsFile.is_open()){
    std::cerr<< "problem opening groups file, giving up!\n";
    exit(1);
  }

  for(int iGroup = 0; iGroup < nGroups; iGroup++){
    for(int iVertex = iGroupStarts[iGroup]; iVertex < iGroupStarts[iGroup+1]; iVertex++){
      groupsFile << iVertex << " " << iGroup << "\n";
    }
  }
  groupsFile.close();
  

  return 0;
}
