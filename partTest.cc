#include "iostream"
#include <string.h>
#include <ctype.h>
#include <vector>
#include <ctime>
#include <cmath>
#include <iomanip>
#include "Graph.hpp"
#include "GraphFactory.hpp"
#include "Sparse.hpp"
#include "GraphPartitioning.hpp"
#include "Runlog.hpp"
#include "Sparse.hpp"
#include "RngWrapper.hpp"

int main(int argc, char **argv){
  bool allOk = true;
	
  std::clock_t  mainStart, localStart;
  
  mainStart = std::clock();
  std::ifstream infile;

  std::string inputDirectory;
  std::string outputDirectory;
  std::string graphName;
  std::string graphNameAndFormat;
  int graphFormat = 0;
  std::string nSplitsString;
  int nSplits = 1;

  infile.open(argv[1]);
  if(!infile.is_open()){
    allOk = false;
    std::cerr << "Couldn't open params file " << argv[1] << "\n";
  }else{
    getline (infile,inputDirectory);
    getline (infile,outputDirectory);
    getline (infile, graphNameAndFormat);
    getline (infile, nSplitsString);
    
    std::stringstream ss(graphNameAndFormat);
    std::string param;
    
    std::getline(ss, param, ' ');
    graphName = param;
      
    std::getline(ss, param, ' ');
    graphFormat = atoi(param.c_str());
    
    nSplits = atoi(nSplitsString.c_str());

    std::cout << "inputDirectory " << inputDirectory << "\n";
    std::cout << "outputDirectory " << outputDirectory << "\n";
    std::cout << "graphName " << graphName << "\n";
    
    std::cout << "Number of partitions: 2^" << nSplits << std::endl;
  }

  char message[200];
  RngWrapper rng;
  GraphFactory graphFactory;
  AdjSMatGraph *mySMATGraph = NULL;
  Sparse<double> laplacian;
  Runlog runlog(outputDirectory,std::string("log.txt"));
  if(allOk){

    std::vector<long> rngSeeds;
    long rngSeed = 1234542;
    
    rng.setSeed(rngSeed);
    
    sprintf (message, "Setting local RNG seed to %lu, first random is %f ", rngSeed, rng.getUniform01Random());
    runlog.writelog(message);
     
    mySMATGraph = (AdjSMatGraph*)graphFactory(2);
    std::string graphFile = std::string(inputDirectory+graphName+std::string(".txt"));

    sprintf (message, "Reading in graph from file %s with format code %d", graphFile.c_str(),graphFormat);
    runlog.writelog(message);

    if(!mySMATGraph-> readFromFile(graphFile.c_str(),graphFormat)){
      std::cerr << "Problem reading data from file " << graphFile << std::endl;
      allOk = false;
    }
    if(allOk){
      std::cout << "Graph is read in at: " << (std::clock() - mainStart) / (double)(CLOCKS_PER_SEC) << " seconds" << std::endl;
    }else{
      std::cout << "Did not succeed in reading graph\n";
    }
  }else{
    sprintf (message, "Partitioning failed");
    runlog.writelog(message);
  }
  
  if(allOk){
    //std::string check("/Users/Martin/Google Drive/HPC/Graphs/partlogs/check.txt");
    //mySMATGraph->writeEdgeListToFile(check);
    
    
    std::cout << "Graph has " << mySMATGraph->nVertices() << " vertices and " << mySMATGraph->nEdges() << " edges\n";
    std::cout << "Density: " << (double)mySMATGraph->nEdges()/(mySMATGraph->nVertices()*mySMATGraph->nVertices()) << "\n";
    
    std::cout << "Splitting into " << pow(2.0,nSplits) << " subdivisions\n";
    
    int nComponents;
    
    std::vector<int> components;
    localStart = std::clock();
    nComponents = mySMATGraph->connectedComps(components);
    std::cout << "Checked connectedness in " << (std::clock() - localStart) / (double)(CLOCKS_PER_SEC) << " seconds" << std::endl;
    
    if(nComponents == 1){
      std::cout << "There is 1 connected component" << std::endl;
    }else{
      std::cout << "There are " << nComponents << " connected components , stopping\n";
      allOk = false;
      return 1;
    }
    
    localStart = std::clock();
    laplacian = mySMATGraph->getLaplacian();
     
    std::vector<int> partition(laplacian.nRows());
    
    GraphPartitioning graphPartitioner(&runlog);
    
    graphPartitioner.partitionGraph(laplacian, partition, nSplits,1, rng);
    
    std::ofstream graphPart;
    std::string pathAndFileName(outputDirectory+graphName+ "_part"+ runlog.to_string(nSplits) + ".txt");
    graphPart.open(pathAndFileName.c_str());
  
    std::cout << "Writing partition to " << outputDirectory << graphName << "_part"+ runlog.to_string(nSplits) + ".txt\n";

    for(int i = 0; i < (int)partition.size(); i++){
	graphPart  << partition [i] << "\n";
    }
    graphPart.close();
 
    delete mySMATGraph;
    
  }
  sprintf (message,  "Time of main: %f seconds",(std::clock() - mainStart) / (double)(CLOCKS_PER_SEC));
  runlog.writelog(message);
  
  std::cout << "Time of main: " << (std::clock() - mainStart) / (double)(CLOCKS_PER_SEC) << " seconds" << std::endl;
  sprintf (message, "***** Done *****");
  runlog.writelog(message);
  
  
  return 0;
}


