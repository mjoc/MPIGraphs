#include "Ising.hpp"
#include "Graph.hpp"
#include "RngWrapper.hpp"
#include "Runlog.hpp"
#include <set>
#include <cmath>
#include <queue>
#include <iostream>

#define OUT_DIR "/home/users/mschpc/2013/oconnm28/Courses/GraphPartitioning/Code/IsingOutput/"

Ising::Ising(AdjSMatGraph* graph, Runlog* runlog, std::string name){
  _graph = graph;
  _colours.assign(graph->nVertices(),UNCOLOURED);
  _states.assign(graph->nVertices(),0);
  _bipartite = untested;
  _runlog = runlog;
  _name = name;
}

void Ising::printColours(){
  
  for(std::vector<int>::const_iterator it = _colours.begin(); it != _colours.end(); ++it){
      _runlog->writelog(_runlog->to_string(*it));
    }
}

bool Ising::isBipartite(int iStartVertex, int iStartColour){
    bool answer = true;
    
    if(_bipartite == untested){
      int iCurrentVertex = iStartVertex;
      
      _colours[iCurrentVertex] = iStartColour;
      
      if(iStartColour==GREEN){
        _greens.insert(iCurrentVertex);
      }else{
        _reds.insert(iCurrentVertex);
      }
      
      std::queue<int> visiting;
      visiting.push(iCurrentVertex);
      
      while(!visiting.empty()){
        iCurrentVertex = visiting.front();
        visiting.pop();
        for(int i = _graph->_rowIndex[iCurrentVertex];i < _graph->_rowIndex[iCurrentVertex+1]; i++){
          if(_colours[_graph->_cols[i]] == UNCOLOURED){
            if(_colours[iCurrentVertex] == GREEN){
              _colours[_graph->_cols[i]] = RED;

              _reds.insert(_graph->_cols[i]);
            }else{
              _colours[_graph->_cols[i]] = GREEN;

              _greens.insert(_graph->_cols[i]);
            }
            visiting.push(_graph->_cols[i]);
          }else{
            if(_colours[_graph->_cols[i]]== _colours[iCurrentVertex]){

              answer = false;
              break;
            }
          }
        }
      }
      
      if(!answer){
        _greens.clear();
        _reds.clear();
      }
    }else{
      if(_bipartite == no){
        answer = false;
      }
      
    }

    return answer;
  }
  
const std::set<int>& Ising::greens(){
    return _greens;
  }
  
  const std::set<int>& Ising::reds(){
    return _reds;
  }
  
  int Ising::nGreen(){
    return (int)_greens.size();
  }

  int Ising::nRed(){
    return (int)_reds.size();
  }
   
  bool Ising::doIsingSimulation(double JoverKBT, int nBurnIn, int nSamples, int stepSize, RngWrapper& rng){
    bool allOk = true;
    int energyDelta = 0;
    double criticalValue;
    std::set<std::set<int> > colours;
    colours.insert(_greens);
    colours.insert(_reds);
    std::vector<double> random01numbers;
    std::cout << "Doing Ising burn-in of " << nBurnIn << std::endl;
    
    int flips = 0;
    
    if(((int)_greens.size() + (int)_reds.size()) != (int)_graph->_nRows){
      
      std::cout << "Problem! green: " << _greens.size() << "  red: " << _reds.size() << " vertices: " << _graph->_nRows <<std::endl;
      allOk = false;
    }else{
      std::vector<int>::iterator it;
      for(it = _states.begin(); it < _states.end(); ++it){
        *it = 1;
      }
      
      double averageFlips = 0;
      for(int i = 0; i < nBurnIn; i++){
        if(i%100 == 0){
          std::cout << "burn in: " << i << " of " << nBurnIn << std::endl;
        }
        random01numbers = rng.getUniform01Randoms(_graph->_nRows);
        flips = 0;
        for(std::set<std::set<int> >::iterator colourIter = colours.begin(); colourIter != colours.end(); ++colourIter){
          for(std::set<int>::iterator it = (*colourIter).begin(); it != (*colourIter).end(); ++it){
            energyDelta = 0;
            for(int j = _graph->_rowIndex[*it];j <  _graph->_rowIndex[*it+1]; j++){
              energyDelta += _states[_graph->_cols[j]];
            }
            energyDelta *= _states[*it];
            criticalValue = exp(-2*JoverKBT*energyDelta);
            
            if(random01numbers[*it] < criticalValue ) {
              _states[*it] *= -1;
              flips++;
            }
          }
        }
        averageFlips += flips;
      }
      
      averageFlips /= nBurnIn;
      std::cout << "Average flips: " << averageFlips << " on " << _graph->_nRows << " vertices (" << averageFlips/_graph->_nRows << "%)"<<std::endl;
      
      
      //File stuff
      std::ofstream samplesFile, rFile, statesStartFile, statesEndFile, paramsFile, flipsFile;
      
      std::string directory = std::string(OUT_DIR);
      std::string samplesFileName(directory + "sample.txt");
      std::string rFileName(directory + "r.txt");
      std::string statesStartFileName(directory + "latticeStart.txt");
      std::string statesEndFileName(directory + "latticeEnd.txt");
      std::string paramsFileName(directory + "params.txt");
      std::string flipsFileName(directory + "flips.txt");
      
      samplesFile.open (samplesFileName.c_str());
      
      rFile.open (rFileName.c_str());
      statesStartFile.open(statesStartFileName.c_str());
      statesEndFile.open(statesEndFileName.c_str());
      paramsFile.open(paramsFileName.c_str());
      flipsFile.open(flipsFileName.c_str());
      
      
      for(std::vector<int>::iterator it = _states.begin(); it != _states.end(); ++it){
        statesStartFile << *it << std::endl;
      }

      std::cout << "Main chain " << std::endl;
      //Main Chain
      int iSample = 0;
      int chainIndex = 0;
      int collectSample = 0;
      int averageSpin = 0;
      std::vector<double> averageSpins(nSamples);
      std::vector<double> flipRate(nSamples);
      
      while(iSample< nSamples){
        if(chainIndex % stepSize == 0)
        {
          //Collect sample
          collectSample = 1;
          averageSpin = 0;
          flipRate[iSample] = (double)flips/(_graph->_nRows);
        }
        flips = 0;
        random01numbers = rng.getUniform01Randoms(_graph->_nRows);
        
        //Do the sampling if nessesary, and always reshuffle
        for(std::set<std::set<int> >::iterator colourIter = colours.begin(); colourIter != colours.end(); ++colourIter){
          
          for(std::set<int>::iterator it = (*colourIter).begin(); it != (*colourIter).end(); ++it){
            //Collect the data before resampling
            if(collectSample == 1){
              averageSpin += _states[*it];
            }
            
            energyDelta = 0;
            for(int j = _graph->_rowIndex[*it];j <  _graph->_rowIndex[*it+1]; j++){
              energyDelta += _states[_graph->_cols[j]];
              
            }
            energyDelta *= _states[*it];
            
            criticalValue = exp(-2*JoverKBT*energyDelta);
            
            if(random01numbers[*it] < criticalValue ) {
              _states[*it] *= -1;
              flips++;
            }
            
          }
        }
        //Aggregate the Sample, if taken
        if(collectSample == 1){
          averageSpins[iSample] = (double)averageSpin/(_graph->_nRows);
          samplesFile << averageSpins[iSample] << std::endl;
          iSample++;
          collectSample = 0;
          
        }
        chainIndex++;
      }
      
      for(std::vector<int>::iterator it = _states.begin(); it != _states.end(); ++it){
        statesEndFile << *it << std::endl;
      }
      
      // Autocorrelation function for t <= 20
      int t  = 20;
      std::vector<double>r(t+1);
      std::vector<double>correlations(t);
      
      
      double overallMean = 0.0;
      for(int i = 0; i < nSamples; i++){
        overallMean += averageSpins[i];
      }
      overallMean /= (double)(nSamples);
      
      for(int lagStep = 0; lagStep <= t; lagStep++){
        double numerator = 0.0;
        int count = 0;
        for(int i = 0; i < nSamples-lagStep; i++){
          numerator  += (averageSpins[i] - overallMean)*(averageSpins[i+lagStep]- overallMean);
          count++;
        }
        
        r[lagStep] = numerator/(count-1) - (overallMean*overallMean);
      }
      for(int lag = 0; lag <= t; lag++){
        rFile << r[lag] << std::endl;
      }
      
      double tauInt = 0.0;
      for(int lag = 1; lag <= t; lag++){
        tauInt += r[lag]/r[0];
      }
      tauInt += 0.5;
      tauInt *= stepSize;
      
      std::cout << "tauInt: " << tauInt << std::endl;
      
      paramsFile << "Graph Size" << _graph->_nRows << std::endl;
      paramsFile << "Green/red Balance" << _greens.size() << "/" << _reds.size() <<std::endl;
      paramsFile << "BurnIn " << nBurnIn << std::endl;
      paramsFile << "StepSize " << stepSize << std::endl;
      paramsFile << "nSamples " << nSamples << std::endl;
      paramsFile << "JoverKBT " << JoverKBT << std::endl;
      //paramsFile << "seed" << seed << std::endl;
      paramsFile.close();
      
      for(int i = 0; i < nSamples; i++){
        flipsFile << flipRate[i]<< std::endl;
      }
      
      //Tidy up
      samplesFile.close();
      rFile.close();
      statesStartFile.close();
      statesEndFile.close();
      flipsFile.close();
      
    }
    return allOk;
  }
