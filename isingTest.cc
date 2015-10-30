#include "mpi.h"
#include <string.h>
#include <ctype.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include "Graph.hpp"
#include "GraphFactory.hpp"
#include "Sparse.hpp"
#include "Runlog.hpp"
#include "DistribGraph.hpp"
#include "RngWrapper.hpp"
#include "Ising.hpp"
#include "DistribIsing.hpp"

enum PartitionType { SPECIFIED, AS_IS, RANDOMISE};
#define PART_TYPE1 "partition"
#define PART_TYPE2 "as_is"
#define PART_TYPE3 "random"



int main(int argc, char **argv){
  int rank, size, allOk = true;  
  MPI_Init (&argc, &argv);	/* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &size);
	
  if(rank == 0){
    std::cout << "There are " << size << " processes running\n";
  }
  

  std::ifstream infile;
  std::string data;
  std::string inputDirectory;
  std::string outputDirectory;
  std::string graphNameAndFormat;
  std::string graphName;
  int graphFormat = 0;
  std::string partitionType;
  
  std::string isingParams;

  double JoverKBT = 0.44;
  int nSamples = 250;
  int nBurnIn = 1000;
  int nSampleStride = 250;
  
  PartitionType partType = AS_IS;
  
  if(rank == 0){
    infile.open(argv[1]);
    if(!infile.is_open()){
      allOk = false;
      std::cerr << "Couldn't open params file " << argv[1] << "\n";
    }else{
      getline (infile,inputDirectory);
      getline (infile,outputDirectory);
      getline (infile, graphNameAndFormat);
      getline (infile,partitionType);
      getline (infile, isingParams);
      
      std::stringstream ss(graphNameAndFormat);
      std::string param;
      
      std::getline(ss, param, ' ');
      graphName = param;
      
      std::getline(ss, param, ' ');
      graphFormat = atoi(param.c_str());
      
      std::string partitionTypeText1(PART_TYPE1);
      std::string partitionTypeText2(PART_TYPE2);
      std::string partitionTypeText3(PART_TYPE3);
      
      if(partitionType.compare(partitionTypeText1) == 0){
        partType = SPECIFIED;
      }else if(partitionType.compare(partitionTypeText2) == 0){
        partType = AS_IS;
      }else if(partitionType.compare(partitionTypeText3) == 0){
        partType = RANDOMISE;
      }
    
      std::cout << "inputDirectory " << inputDirectory << "\n";
      std::cout << "outputDirectory " << outputDirectory << "\n";
      std::cout << "graphName " << graphName << "\n";
      if(partType == SPECIFIED){
        std::cout << "Partition Type " << PART_TYPE1<<"\n";
      }else if(partType == AS_IS){
        std::cout << "Partition Type " << PART_TYPE2 <<"\n";
      }else if(partType == RANDOMISE){
        std::cout << "Partition Type " << PART_TYPE3 <<"\n";
      }else{
        std::cout << "Partition Type: default\n";
      }
      std::cout << "Ising parameters: " << isingParams << std::endl;
    }
  }
  
  MPI_Bcast(&allOk, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(!allOk ){
    if (rank == 0) {
      std::cerr << "Error: Program terminating\n";
    }
    MPI_Finalize();
    exit(1);
  }
  
  int paramSizes[3];
  if(rank == 0){
    paramSizes[0] = (int)inputDirectory.length()+1;
    paramSizes[1] = (int)outputDirectory.length()+1;
    paramSizes[2] = (int)graphName.length()+1;
  }
  MPI_Bcast(paramSizes, 3, MPI_INT, 0, MPI_COMM_WORLD);
  
  char inDir[paramSizes[0]];
  char outDir[paramSizes[1]];
  char name[paramSizes[2]];
  
  
  if(rank == 0){
    strncpy(inDir,inputDirectory.c_str(),paramSizes[0]);
    strncpy(outDir,outputDirectory.c_str(),paramSizes[1]);
    strncpy(name,graphName.c_str(),paramSizes[2]);
  }
  
  MPI_Bcast(inDir, paramSizes[0], MPI_CHAR, 0, MPI_COMM_WORLD);

  MPI_Bcast(outDir, paramSizes[1], MPI_CHAR, 0, MPI_COMM_WORLD);

  MPI_Bcast(name, paramSizes[2], MPI_CHAR, 0, MPI_COMM_WORLD);

  double isingParameters[4];
  
  if(rank == 0){
    std::stringstream ss(isingParams);
    std::string param;
    
    std::getline(ss, param, ' ');
    JoverKBT = atof(param.c_str());
    
    std::getline(ss, param, ' ');
    nSamples = atof(param.c_str());
    
    std::getline(ss, param, ' ');
    nBurnIn = atoi(param.c_str());
    
    std::getline(ss, param, ' ');
    nSampleStride = atoi(param.c_str());
    
    isingParameters[0] = JoverKBT;
    isingParameters[1] = nSamples;
    isingParameters[2] = nBurnIn;
    isingParameters[3] = nSampleStride;
    std::cout << "nSamples : " <<  nSamples << std::endl;
    std::cout << "nBurnIn :" << nBurnIn << std::endl;
    std::cout << "nSampleStride :" << nSampleStride << std::endl;
    std::cout << "JoverKBT :" << JoverKBT << std::endl;
  }
  MPI_Bcast(isingParameters, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if(rank != 0){
     JoverKBT = isingParameters[0];
     nSamples = isingParameters[1];
     nBurnIn = isingParameters[2];
     nSampleStride = isingParameters[3];
  }
  
  if(rank != 0){
    
    inputDirectory = std::string(inDir);
    outputDirectory = std::string(outDir);
    graphName = std::string(name);
    
  }
	
  char message[200];
  std::string filename(graphName + std::string("_log.txt"));
  Runlog runlog(outputDirectory,filename,rank);
  
  RngWrapper rng;
  rng.setSeed(451239);
  std::vector<long> rngSeeds;
  long rngSeed;
  if(rank == 0){
    rngSeeds = rng.getSeeds(size);
  }
  
  MPI_Scatter(&rngSeeds[0],1,MPI_LONG,&rngSeed,1, MPI_LONG, 0 , MPI_COMM_WORLD);
  
  rng.setSeed(rngSeed);
  
  sprintf (message, "Setting local RNG seed to %lu, first random is %f ", rngSeed, rng.getUniform01Random());
  runlog.writelog(message);
 	
  GraphFactory graphFactory;
  AdjSMatGraph *mySMATGraph = 0;
  Sparse<double> r1laplacian;
  if(rank == 0){

    sprintf (message, "Reading in graph from file %s with format code %d", graphName.c_str(), graphFormat);
    runlog.writelog(message);
    
    mySMATGraph = (AdjSMatGraph *)graphFactory(2);
    std::string graphFile = std::string(inputDirectory+graphName+std::string(".txt"));
    if(!mySMATGraph->readFromFile(graphFile.c_str(),graphFormat)){
      std::cerr << "Problem reading data from file " << graphFile << std::endl;
      allOk = false;
    }
  }
  MPI_Bcast(&allOk, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(!allOk ){
    if (rank == 0) {
      std::cerr << "Error: Program terminating\n";
    }
    MPI_Finalize();
    exit(1);
  }

  bool doSerial = false;
  if(rank == 0 && doSerial){
    std::cout << "Graph has " << mySMATGraph->nVertices() << " vertices and " << mySMATGraph->nEdges() << " edges\n";
    std::clock_t timerStart = std::clock();
    int nComponents;
    std::vector<int> components;
    
    nComponents = mySMATGraph->connectedComps(components);
    
    if(nComponents == 1){
      std::cout << "There is 1 connected component\n";
    }else{
      std::cout << "There are " << nComponents << " connected components\n";
    }
    
    if(nComponents == 1){
      Ising isingModel(mySMATGraph, &runlog, graphName);
      
      if(isingModel.isBipartite(0,1)){
        std::cout << "This graph is bipartite" << std::endl;
        std::cout << "Green: " << isingModel.nGreen() << "  red: " << isingModel.nRed() << std::endl;
        
        isingModel.doIsingSimulation(JoverKBT, nBurnIn, nSamples, nSampleStride, rng);
			}else{
        std::cout << "This graph is not bipartite" << std::endl;
      }
    }
    std::cout << "Serial run takes " << (std::clock()-timerStart)/(double)(CLOCKS_PER_SEC) << " seconds";
  }
  bool doParallel = true;
  if(doParallel){
    sprintf (message, "Distributing Graph");
    runlog.writelog(message);
    
    DistribGraph distGraph(MPI_COMM_WORLD, rank, size, &runlog);
    
    if(rank == 0){
      if(partType == SPECIFIED){
	std::cout << "Partition specified\n";
	std::string partitionFile(inputDirectory+graphName+std::string("_part.txt"));
	distGraph.distributeGraphRankZero(*mySMATGraph, partitionFile.c_str());
      }else if(partType == RANDOMISE){
	int baseSize = (int)mySMATGraph->nVertices()/size;
	int residual = (int)mySMATGraph->nVertices() - size * baseSize;
	std::vector<int> nLocalVerticesByRank(size);
	{
	  int i;
	  for(i = 0; i < size; i++){
	    nLocalVerticesByRank[i] = baseSize;
	  }
	  
	  i = 0;
	  while(residual > 0 && i  < size){
	    nLocalVerticesByRank[i]++;
	    i++;
	    residual--;
	  }
	}
	{
	  std::vector<int> randomOrder = rng.getRandomPermutation(mySMATGraph->nVertices());
	  std::vector<int> randomRanks(mySMATGraph->nVertices());
	  int iVertex = 0;
	  for(int iRank = 0; iRank < size; iRank++){
	    for(int iVertexInRank = 0; iVertexInRank < nLocalVerticesByRank[iRank]; iVertexInRank++){
	      randomRanks[randomOrder[iVertex]] = iRank;
	      iVertex++;
	    }
	  }
	  distGraph.distributeGraphRankZero(*mySMATGraph, randomRanks);
	}
      }else{
	distGraph.distributeGraphRankZero(*mySMATGraph);
      }
    
    }else{
      distGraph.receiveGraphRankNonZero();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    sprintf (message, "Graph distributed");
    runlog.writelog(message);
    
    int externalEdges = distGraph.nExternalEdges();
    int interProcessEdges = 0;
    
    MPI_Reduce(&externalEdges, &interProcessEdges, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if(rank == 0){
      sprintf(message, "External Edges: %d of %d total (%0.2f)",interProcessEdges,mySMATGraph->nEdges(),(double)interProcessEdges/mySMATGraph->nEdges());
      runlog.writelog(message);
    }
    
    sprintf (message, "Setting up distributed ising model");
    runlog.writelog(message);
    DistribIsing distIsingModel(rank, size, MPI_COMM_WORLD, &distGraph, &runlog,graphName, outputDirectory);
    
    std::clock_t  mainStart = 0, localStart = 0;
    if(rank == 0){
      mainStart = std::clock();
    }
    sprintf(message, "Checking Bipartite: ");
    runlog.writelog(message);
    localStart = std::clock();
    bool isBipart = distIsingModel.isBipartite(0);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0){
      sprintf (message, "Time of Bipartite: %f seconds", (std::clock() - localStart) / (double)(CLOCKS_PER_SEC));
      
      runlog.writelog(message);
    }
    if(isBipart){
      sprintf(message, "Bipartite checked and is ok, colours set up\n doing running Ising model");
      runlog.writelog(message);
      localStart = std::clock();
      distIsingModel.doIsingSimulation(JoverKBT, nBurnIn, nSamples, nSampleStride, rng);
      MPI_Barrier(MPI_COMM_WORLD);
      if(rank == 0){
	sprintf (message,"Time of Ising: %f seconds",  (std::clock() - localStart) / (double)(CLOCKS_PER_SEC));
	runlog.writelog(message);
      }
    }else{
      sprintf(message, "Bipartite checked and NOT ok!");
      runlog.writelog(message);
    }
    if(rank == 0){
      sprintf (message,"Time of algorithm: %f seconds",  (std::clock() - mainStart) / (double)(CLOCKS_PER_SEC));
      runlog.writelog(message);
      
    }
    
    if(rank == 0){
      delete mySMATGraph;
    }
  }  
  sprintf (message, "***** Done: Rank:  %d", rank);
  runlog.writelog(message);
  
  MPI_Finalize();
  std::cout << "Rank "<< rank << " is done!\n";
  
  return 0;
}



