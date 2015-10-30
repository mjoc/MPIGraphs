#include "DistribIsing.hpp"
#include "Runlog.hpp"
#include "RngWrapper.hpp"
#include <iostream>


void DistribIsing::statusMessage(char* message){
  if(_runlog != NULL){
    _runlog->writelog(message);
  }else{
    std::cout << message << std::endl;
  }
}

void DistribIsing::sendAndReceive(std::vector<int> sendVector){
  for(int i = 0; i < (int)_sendData.size(); i++){
    _sendData[i] = sendVector[_sendDataInVertex[i]];
  }
  
  MPI_Alltoallv(&_sendData[0],&_sendDataNByRank[0],&_sendDataRankDisplacements[0], MPI_INT,&_recvData[0],&_recvDataNByRank[0],&_recvDataRankDisplacements[0], MPI_INT, _communicator);

}


void DistribIsing::sendAndReceive(std::vector<int> sendVector, int colour){
  std::clock_t timerStart;
  std::clock_t timerPackUnpackStart;
  std::clock_t timerAll2All;

  timerStart = std::clock();
  if(colour == GREEN){
    
    _dataVolumeOut += (double)_sendDataGreen.size();

    timerPackUnpackStart = std::clock();
    for(int i = 0; i < (int)_sendDataGreen.size(); i++){
      _sendDataGreen[i] = sendVector[_sendDataGreenInVertex[i]];
    }
    _timerPackUnpack += (std::clock()-timerPackUnpackStart)/(double)(CLOCKS_PER_SEC);

    timerAll2All = std::clock();
    MPI_Alltoallv(&_sendDataGreen[0],&_sendDataGreenNByRank[0],&_sendDataGreenRankDisplacements[0], MPI_INT,&_recvDataRed[0],&_recvDataRedNByRank[0],&_recvDataRedRankDisplacements[0], MPI_INT, _communicator);
    _timerAll2All = (std::clock()-timerAll2All)/(double)(CLOCKS_PER_SEC);

    _dataVolumeIn += (double)_recvDataRed.size();

    timerPackUnpackStart = std::clock();
    for(int i = 0; i < (int)_recvDataRed.size(); i++){
      _recvData[_recvDataRedToFullIndex[i]] = _recvDataRed[i];
     }
    _timerPackUnpack += (std::clock()-timerPackUnpackStart)/(double)(CLOCKS_PER_SEC);

  }
  
  if(colour == RED){
    _dataVolumeOut += (double)_sendDataRed.size();

    timerPackUnpackStart = std::clock();
    for(int i = 0; i < (int)_sendDataRed.size(); i++){
      _sendDataRed[i] = sendVector[_sendDataRedInVertex[i]];
    }
    _timerPackUnpack += (std::clock()-timerPackUnpackStart)/(double)(CLOCKS_PER_SEC);
     
    timerAll2All = std::clock();
    MPI_Alltoallv(&_sendDataRed[0],&_sendDataRedNByRank[0],&_sendDataRedRankDisplacements[0], MPI_INT,&_recvDataGreen[0],&_recvDataGreenNByRank[0],&_recvDataGreenRankDisplacements[0], MPI_INT, _communicator);
    _timerAll2All = (std::clock()-timerAll2All)/(double)(CLOCKS_PER_SEC);

    _dataVolumeIn += (double)_recvDataGreen.size();

    timerPackUnpackStart = std::clock();
    for(int i = 0; i < (int)_recvDataGreen.size(); i++){
      _recvData[_recvDataGreenToFullIndex[i]] = _recvDataGreen[i];
    }
    _timerPackUnpack += (std::clock()-timerPackUnpackStart)/(double)(CLOCKS_PER_SEC);
     
  }
  _timer += (std::clock()-timerStart)/(double)(CLOCKS_PER_SEC);
}



////////////////////////////////////////////////
bool DistribIsing::isBipartite(int iStartVertex){
  char message[200];
  bool allOk = true;
  bool isBipart = true;

  if(_mpi_rank==0){
    int nComponents;
    std::vector<int> components;
    
    nComponents = _graph->connectedComps(components);
    
    std::vector<int> componentCounts(nComponents);
    std::vector<int> componentColours(nComponents);
    
    for(std::vector<int>::const_iterator it = components.begin(); it != components.end(); ++it){
      componentCounts[*it]++;
    }
    
    int maxComponentSize = 0;
    int biggestLocalComponent = -1;
    {
      int i = 0;
      for(std::vector<int>::const_iterator it = componentCounts.begin(); it != componentCounts.end(); ++it){
        if (*it > maxComponentSize){
          maxComponentSize = *it;
          biggestLocalComponent = i;
        }
        i++;
      }
    }
    
    {
      int iVertexInBiggestComponent = 0;
      while(components[iVertexInBiggestComponent] != biggestLocalComponent){
        iVertexInBiggestComponent++;
      }
      isBipart = Ising::isBipartite(iVertexInBiggestComponent, GREEN);
    }
    
    std::vector<int> biggestComponents(_mpi_size);
    
    
  }
  bool keepGoing = true, notBipartite = false;
  
  int localFinishedFlag = 0;
  int globalFinishedCount = 0;
  while(keepGoing){
    sendAndReceive(_colours);
    for(int i = 0; i < (int)_recvData.size(); i++){
      for(int j = 0; j < (int)_recvDataInVertices[i].size(); j++){
        int iReceivingVertex =  _recvDataInVertices[i][j];
        if(_recvData[i] != 0){
          int newColour = _recvData[i]==GREEN?RED:GREEN;
          
          if(_colours[iReceivingVertex] == 0 ){
            
            int iColour = _recvData[i]==GREEN?RED:GREEN;
            isBipart = Ising::isBipartite(iReceivingVertex, iColour);
            if(!isBipart){
              notBipartite = true;
              sprintf(message, "Component with vertex %d is not bipartite", _distribGraph->_globalVertices[ iReceivingVertex]);
              _runlog->writelog(message);
            }
          }else if (_colours[iReceivingVertex] != newColour ){
            notBipartite = true;
          }
        }
      }
    }
    
    localFinishedFlag = 1;
    if(!notBipartite){
      for(int i = 0; i < (int)_colours.size(); i++){
        if(_colours[i]==UNCOLOURED){
          localFinishedFlag = 0;
          break;
        }
      }
    }else{
      localFinishedFlag = -(_mpi_size + 1);
    }
    
    MPI_Allreduce(&localFinishedFlag,&globalFinishedCount,1,MPI_INT,MPI_SUM,_communicator);
    
    if(globalFinishedCount == _mpi_size || globalFinishedCount < 0){
      if(globalFinishedCount < 0){
        isBipart = false;
      }
      keepGoing = false;
      if(!isBipart){
        sprintf(message, "Doesn't seem to be bipartite, localFinishedFlag is %d", localFinishedFlag);
        _runlog->writelog(message);
      }

    }

  }
  
  int nRed = 0, nGreen = 0, nGlobalGreen, nGlobalRed;
  for(std::vector<int>::const_iterator it = _colours.begin(); it != _colours.end(); it++){
    if(*it == GREEN){nGreen++;}else if(*it == RED){nRed++;}
  }
  MPI_Reduce(&nGreen, &nGlobalGreen, 1, MPI_INT, MPI_SUM, 0, _communicator);
  MPI_Reduce(&nRed, &nGlobalRed, 1, MPI_INT, MPI_SUM, 0, _communicator);
  if(_mpi_rank==0){
    sprintf(message, "There are %d green and %d red vertices, total %d", nGlobalGreen, nGlobalRed, _distribGraph->nGlobalVertices());
    _runlog->writelog(message);
    
  }
  
  if(isBipart){
    _sendDataGreenNByRank.resize(_mpi_size);
    _sendDataRedNByRank.resize(_mpi_size);
    _sendDataGreenRankDisplacements.resize(_mpi_size);
    _sendDataRedRankDisplacements.resize(_mpi_size);
   
    int iVertex = 0;
    int nSendDataGreen = 0, nSendDataRed = 0;
    for(int iRank = 0; iRank < _mpi_size; iRank++){
      for(int iVertexInRank = 0; iVertexInRank < _sendDataNByRank[iRank]; iVertexInRank++){
        
        if(_colours[_sendDataInVertex[iVertex]]==GREEN){
          _sendDataGreenInVertex.push_back(_sendDataInVertex[iVertex]);
          _sendDataGreenInVertexGlobal.push_back(_sendDataInVertexGlobal[iVertex]);
          _sendDataGreenNByRank[iRank]++;
        }else{
          _sendDataRedInVertex.push_back(_sendDataInVertex[iVertex]);
          _sendDataRedInVertexGlobal.push_back(_sendDataInVertexGlobal[iVertex]);
          _sendDataRedNByRank[iRank]++;
        }
        iVertex++;
      }
      if(iRank == 0){
        _sendDataGreenRankDisplacements[0] = 0;
      }else{
        _sendDataGreenRankDisplacements[iRank] = _sendDataGreenRankDisplacements[iRank-1]+_sendDataGreenNByRank[iRank-1];
      }
      if(iRank == 0){
        _sendDataRedRankDisplacements[0] = 0;
      }else{
        _sendDataRedRankDisplacements[iRank] = _sendDataRedRankDisplacements[iRank-1]+_sendDataRedNByRank[iRank-1];
      }
      nSendDataGreen +=  _sendDataGreenNByRank[iRank];
      nSendDataRed += _sendDataRedNByRank[iRank];
    }
    
    _sendDataGreen.resize(nSendDataGreen);
    _sendDataRed.resize(nSendDataRed);
    
    _recvDataGreenNByRank.resize(_mpi_size);
    _recvDataRedNByRank.resize(_mpi_size);
    _recvDataGreenRankDisplacements.resize(_mpi_size);
    _recvDataRedRankDisplacements.resize(_mpi_size);
    
    iVertex = 0;
    int nRecvDataGreen = 0, nRecvDataRed = 0;
    for(int iRank = 0; iRank < _mpi_size; iRank++){
      for(int iVertexInRank = 0; iVertexInRank < _recvDataNByRank[iRank]; iVertexInRank++){
        if(_colours[_recvDataInVertices[iVertex][0]]==GREEN){
          _recvDataGreenToFullIndex.push_back(iVertex);
          _recvDataGreenNByRank[iRank]++;
        }else{
          _recvDataRedToFullIndex.push_back(iVertex);
          _recvDataRedNByRank[iRank]++;
        }
        
        iVertex++;
      }
      if(iRank == 0){
        _recvDataGreenRankDisplacements[0] = 0;
        _recvDataRedRankDisplacements[0] = 0;
      }else{
        _recvDataGreenRankDisplacements[iRank] =  _recvDataGreenRankDisplacements[iRank-1] + _recvDataGreenNByRank[iRank-1];;
        _recvDataRedRankDisplacements[iRank] = _recvDataRedRankDisplacements[iRank-1] + _recvDataRedNByRank[iRank-1];
      }
      nRecvDataGreen +=  _recvDataGreenNByRank[iRank];
      nRecvDataRed += _recvDataRedNByRank[iRank];

    }
    _recvDataGreen.resize(nRecvDataGreen);
    _recvDataRed.resize(nRecvDataRed);
    
  }
  
  sprintf(message, "There are %d green and %d red colours on the external edges", nGreen,nRed);
  _runlog->writelog(message);
  
 
  std::vector<int> check(_mpi_size);
  
  
  MPI_Alltoall(&_sendDataGreenNByRank[0],1, MPI_INT,&check[0],1, MPI_INT, _communicator);
  for(int i = 0; i < _mpi_size; i++){
    if(check[i] != _recvDataRedNByRank[i]){
      sprintf(message, "Expecting %d from green, got %d for rank %d", _recvDataRedNByRank[i],check[i], i);
      _runlog->writelog(message);
      allOk = false;
    }
  }
  if(allOk){
    sprintf(message, "For green send, expected number received");
    _runlog->writelog(message);
    
  }
  allOk = true;
  
  MPI_Alltoall(&_sendDataRedNByRank[0],1, MPI_INT,&check[0],1, MPI_INT, _communicator);
  for(int i = 0; i < _mpi_size; i++){
    if(check[i] != _recvDataGreenNByRank[i]){
      sprintf(message, "Expecting %d from red, got %d for rank %d", _recvDataGreenNByRank[i],check[i], i);
      _runlog->writelog(message);
      allOk = false;
    }
  }
  if(allOk){
    sprintf(message, "For red send, expected number received");
    _runlog->writelog(message);
    
  }

  sprintf(message, "# Send & Recv Data by Rank SG RG : SR RR");
    _runlog->writelog(message);
    int totalRed= 0, totalGreen = 0;  
  for(int i = 0; i < _mpi_size; i++){
    sprintf(message, "%d %d  %d : %d %d",i,_sendDataGreenNByRank[i],_recvDataGreenNByRank[i],_sendDataRedNByRank[i],_recvDataRedNByRank[i]);
    totalGreen += _sendDataGreenNByRank[i];
    totalRed+= _sendDataRedNByRank[i];
    _runlog->writelog(message);
  }
  sprintf(message, "# Send Data Total Green and Red: %d %d",totalGreen,totalRed);
  _runlog->writelog(message);

  if(allOk){
    sprintf (message, "Colour Edge communication looks good");
  }else{
    sprintf (message, "************Something wrong Colour Edge communication");
  }
  _runlog->writelog(message);
    
  return isBipart;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool DistribIsing::doIsingSimulation(double JoverKBT, int nBurnIn, int nSamples, int stepSize, RngWrapper& rng){
  char message[200];
  bool allOk = true;
  double energyDelta = 1.0, criticalValue;
  std::map<int,std::set<int> > colours;
  colours.insert(std::pair<int,std::set<int> > (GREEN,_greens));
  colours.insert(std::pair<int,std::set<int> > (RED,_reds));
  std::vector<double> random01numbers;
  _timer = 0.0;
  _timerPackUnpack = 0.0;
  _timerAll2All = 0.0;
  _dataVolumeIn = 0.0;
  _dataVolumeOut = 0.0;

  int flips = 0;
  double averageFlips = 0;
  if(((int)_greens.size() + (int)_reds.size()) != (int)_distribGraph->nVertices()){
    
    sprintf(message, "Problem! green: %ld  red: %ld vertices: %d",  _greens.size(), _reds.size(), _distribGraph->nVertices());
    _runlog->writelog(message);

    std::cout << "Problem! green: " << _greens.size() << "  red: " << _reds.size() << " vertices: " << _distribGraph->nVertices() <<std::endl;
    allOk = false;
  }else{
    for(int i = 0; i < (int)_states.size(); ++i){
      _states[i] = 1;
    }
    

    sprintf(message,  "Checking the communication");
    _runlog->writelog(message);
    
  
    sprintf(message,  "Checking edge communication");
    _runlog->writelog(message);

    std::vector<int> localVerticesWithOffset(_distribGraph->nVertices());
    for(int i = 0 ; i < (int)localVerticesWithOffset.size(); i++){
      localVerticesWithOffset[i] = i + _distribGraph->vertexOffset();
    }
    
    sendAndReceive(localVerticesWithOffset,GREEN);
    sendAndReceive(localVerticesWithOffset,RED);
    
    for(std::map<int,std::set<int> >::iterator colourIter = colours.begin(); colourIter != colours.end(); ++colourIter){
      for(std::set<int>::iterator sameColourIter = (colourIter->second).begin(); sameColourIter != (colourIter->second).end(); ++sameColourIter){
        
        // Compile required communications
        std::set<int> connections;
        for(int j = _distribGraph->extendedAdj._rowIndex[*sameColourIter]; j < _distribGraph->extendedAdj._rowIndex[*sameColourIter+1]; j++){
          connections.insert(_distribGraph->extendedAdj._cols[j]);
        }
        
        // Internal communications
        for(int j = _distribGraph->_rowIndex[*sameColourIter]; j  < _distribGraph->_rowIndex[(*sameColourIter)+1]; j++){
          if(connections.size() > 0){
            connections.erase(_distribGraph->vertexOffset() +_distribGraph->_cols[j]);
          }else{
            allOk = false;
            sprintf(message,  "***There is a mismatch with expected vertex communication");
            _runlog->writelog(message);
          }
        }
        
        // External communications
        for(int j = 0; j < (int)_recvDataIndexByInVertex[*sameColourIter].size();j++){
          if(connections.size() > 0){
            connections.erase( _recvData[_recvDataIndexByInVertex[*sameColourIter][j]]);
          }else{
            allOk = false;
            sprintf(message,  "***There is a mismatch with expected vertex communication");
            _runlog->writelog(message);
          }
          
        }
        if(connections.size() > 0){
          allOk = false;
          sprintf(message,  "***Not all expected vertex communications received");
          _runlog->writelog(message);
        }
      }
      if(!allOk){break;}
    }
  }
  if(allOk){
    sprintf(message,  "Edge communication seems ok");
    _runlog->writelog(message);
    
  }
  
  if(allOk){
    sprintf(message,  "Checking colour match");
    
    sendAndReceive(_colours,GREEN);
    sendAndReceive(_colours,RED);
    
    for(std::map<int,std::set<int> >::iterator colourIter = colours.begin(); colourIter != colours.end(); ++colourIter){
      for(std::set<int>::iterator sameColourIter = (colourIter->second).begin(); sameColourIter != (colourIter->second).end(); ++sameColourIter){
        for(int j = _distribGraph->_rowIndex[*sameColourIter];j <  _distribGraph->_rowIndex[*sameColourIter+1]; j++){
          if(_colours[_distribGraph->_cols[j]] == colourIter->first){
            sprintf(message,  "Internal vertex is not of expected colour %d",colourIter->first);
            _runlog->writelog(message);
            allOk  = false;
            break;
          }
          for(int j = 0; j < (int)_recvDataIndexByInVertex[*sameColourIter].size();j++){
            if(_recvData[_recvDataIndexByInVertex[*sameColourIter][j]] == colourIter->first){
              sprintf(message,  "External vertex is not of the expected colour");
              _runlog->writelog(message);
              allOk  = false;
              break;
            }
          }
          if(!allOk){break;}
        }
        if(!allOk){break;}
      }
      if(!allOk){break;}
    }
    if(allOk){
      sprintf(message,  "Colour communication seems ok");
      _runlog->writelog(message);
      
    }
    
    
    if(allOk){
      sprintf(message,  "Doing burning of %d",nBurnIn);
      _runlog->writelog(message);
      
      
      for(int i = 0; i < nBurnIn; i++){
        random01numbers = rng.getUniform01Randoms(_distribGraph->nVertices());
        flips = 0;
        for(std::map<int,std::set<int> >::iterator colourIter = colours.begin(); colourIter != colours.end(); ++colourIter){
          for(std::set<int>::iterator sameColourIter = (colourIter->second).begin(); sameColourIter != (colourIter->second).end(); ++sameColourIter){
            energyDelta = 0;
            
            for(int j = _distribGraph->_rowIndex[*sameColourIter];j <  _distribGraph->_rowIndex[*sameColourIter+1]; j++){
              energyDelta += _states[_distribGraph->_cols[j]];
            }
            for(int j = 0; j < (int)_recvDataIndexByInVertex[*sameColourIter].size();j++){
              energyDelta += _recvData[_recvDataIndexByInVertex[*sameColourIter][j]];
            }
            energyDelta *= _states[*sameColourIter];
            
            criticalValue = exp(-2*JoverKBT*energyDelta);
            
            if(random01numbers[*sameColourIter] < criticalValue ) {
              _states[*sameColourIter] *= -1;
              flips++;
            }
            
          }
          sendAndReceive(_states,colourIter->first);
        }
        averageFlips += flips;
      }
      
      averageFlips /= nBurnIn;
  
      std::ofstream samplesFile;
      char fileName[200];
      
      snprintf(fileName,200,"r%d_%s_sample.txt", _mpi_rank,_name.c_str());
      std::string fullPathAndFile(_directory + std::string(fileName));
      samplesFile.open(fullPathAndFile.c_str());
      
 
      
      //Main Chain
      int iSample = 0;
      int chainIndex = 0;
      int collectSample = 0;
      int aggregateSpin = 0;
      //std::vector<double> averageSpins(nSamples);
      std::vector<double> flipRate(nSamples);
      
      sprintf(message, "Doing Main Chain to get %d samples with stride %d ", nSamples,stepSize);
      _runlog->writelog(message);
      
      while(iSample< nSamples){
        if(chainIndex % stepSize == 0)
        {
          //Collect sample
          collectSample = 1;
          aggregateSpin = 0;
          flipRate[iSample] = (double)flips/(_distribGraph->nVertices());
          

        }
        flips = 0;
        random01numbers = rng.getUniform01Randoms(_distribGraph->nVertices());
        
        for(std::map<int,std::set<int> >::iterator colourIter = colours.begin(); colourIter != colours.end(); ++colourIter){
          
          for(std::set<int>::iterator sameColourIter = (colourIter->second).begin(); sameColourIter != (colourIter->second).end(); ++sameColourIter){
            
            
            if(collectSample == 1){
              aggregateSpin += _states[*sameColourIter];
            }
            energyDelta = 0.0;
            for(int j = _distribGraph->_rowIndex[*sameColourIter];j <  _distribGraph->_rowIndex[*sameColourIter+1]; j++){
              energyDelta += _states[_distribGraph->_cols[j]];
            }
            for(int j = 0; j < (int)_recvDataIndexByInVertex[*sameColourIter].size();j++){
              energyDelta += _recvData[_recvDataIndexByInVertex[*sameColourIter][j]];
            }
            energyDelta *= _states[*sameColourIter];
            
            criticalValue = exp(-2*JoverKBT*energyDelta);
            
            if(random01numbers[*sameColourIter] < criticalValue ) {
              _states[*sameColourIter] *= -1;
              flips++;
            }
            
          }
          sendAndReceive(_states,colourIter->first);
        }
        
        //Aggregate the Sample, if taken
        if(collectSample == 1){
          samplesFile << aggregateSpin << std::endl;
          iSample++;
          collectSample = 0;
          
        }
        chainIndex++;
      }
      
      //Tidy up
      samplesFile.close();

    }
    
  }
  sprintf(message, "Total Send and Receives time: %f seconds", _timer);
  _runlog->writelog(message);
  sprintf(message, "Pack and unpack time: %f seconds", _timerPackUnpack);
  _runlog->writelog(message);
   sprintf(message, "All2All time: %f seconds", _timerAll2All);
  _runlog->writelog(message);

  sprintf(message, "Data out: %f ", _dataVolumeOut) ;
  _runlog->writelog(message);
  sprintf(message, "Data in: %f ", _dataVolumeIn) ;
  _runlog->writelog(message);
  double totalDataOut = 0.0, totalDataIn = 0.0;
  MPI_Reduce(&_dataVolumeOut,&totalDataOut,1,MPI_DOUBLE,MPI_SUM,0,_communicator);
  if(_mpi_rank==0){
    sprintf(message, "Total Sends across ranks: %d (ints)", (int)totalDataOut);
     _runlog->writelog(message);
  }
  MPI_Reduce(&_dataVolumeIn,&totalDataIn,1,MPI_DOUBLE,MPI_SUM,0,_communicator);
  if(_mpi_rank==0){
    sprintf(message, "Total Recvs across ranks: %d (ints)", (int)totalDataIn);
     _runlog->writelog(message);
  }

  return allOk;
}

void DistribIsing::printReceiveDataInternalVertices(){
  char message[200];
  snprintf (message,200, "Where the external data is arriving at");
  _runlog->writelog(message);
  
  for(int i = 0; i < (int)_recvDataInVertices.size(); i++){
    snprintf (message,200, "Vertex %d: ",_sendDataExVertexGlobal[i]);
    for(std::vector<int>::iterator it = _recvDataInVertices[i].begin(); it != _recvDataInVertices[i].end();++it){
      snprintf (message,200, "%s %d",message,(int)*it);
    }
    _runlog->writelog(message);
  }
}
