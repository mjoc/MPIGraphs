#include "mpi.h"
#include <vector>
#include <iostream>
#include <ctime>
#include "RngWrapper.hpp"
#include <unistd.h>

int main(int argc, char **argv){
  int rank, size;
  MPI_Init (&argc, &argv);	/* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  int total = 64;
  int N = 250000;

  unsigned int microseconds;

  std::vector<int> sendData(total);
   std::vector<int> recvData(total);
  std::vector<int> sendDataNByRank(size);
  std::vector<int> recvDataNByRank(size);
  std::vector<int> sendDataRankDisplacements(size);
  std::vector<int> recvDataRankDisplacements(size);
  int totalVolume = 0;
  sendDataRankDisplacements[0] = 0;
  recvDataRankDisplacements[0] = 0;
  for(int i = 1; i <= 2; i++){
    sendDataNByRank[(rank+(7*i))%size] = 32; 
    recvDataNByRank[(rank+size-(7*i))%size] = 32; 
    totalVolume +=sendDataNByRank[(rank+(7*i))%size];
  }
  for(int i = 1; i < size; i++){
    sendDataRankDisplacements[i] = sendDataRankDisplacements[i-1] + sendDataNByRank[i-1];
    recvDataRankDisplacements[i] = recvDataRankDisplacements[i-1] + recvDataNByRank[i-1];
  }

  /*
  if(rank == 0){
    std::cout << "Data per rank send, cum send, recv, cum recv\n";
    for(int i = 0; i < size; i++){
      std::cout << sendDataNByRank[i] << " " << sendDataRankDisplacements[i] << " " << recvDataNByRank[i] << " " << recvDataRankDisplacements[i] << "\n"; 
    }
  }
  */
  std::clock_t timerStart = std::clock();
  for(int i = 0; i < N; i++){
    microseconds = rand() % 50;
    usleep(microseconds);
    
    MPI_Alltoallv(&sendData[0],&sendDataNByRank[0],&sendDataRankDisplacements[0], MPI_INT,&recvData[0],&recvDataNByRank[0],&recvDataRankDisplacements[0], MPI_INT, MPI_COMM_WORLD);
  }
  if(rank==0){
    std::cout << "2 ranks, high volume (" << totalVolume << ") each: " << (std::clock()-timerStart)/(double)(CLOCKS_PER_SEC) << std::endl;
  }
  
  for(int i = 0; i < size; i++){
    sendDataNByRank[i] = 0; 
    recvDataNByRank[i] = 0; 
    sendDataRankDisplacements[i] = 0;
    recvDataRankDisplacements[i] = 0;
  }
  

  totalVolume = 0;
  int ranks[] = {0,2,4,5,7,9,11,13};
  bool activerank = false;
  for(int i = 0; i < 8; i++){
    if(ranks[i]==rank){
      activerank = true;
    }
  }

  if(activerank){
    for(int i = 0; i < 8; i++){
      sendDataNByRank[ranks[i]] = 6; 
      recvDataNByRank[ranks[(8+i-1)%8]] = 6; 
      totalVolume += sendDataNByRank[ranks[i]];
    }
    sendDataNByRank[rank] = 0; 
    recvDataNByRank[rank] = 0; 
    for(int i = 1; i < size; i++){
      sendDataRankDisplacements[i] = sendDataRankDisplacements[i-1] + sendDataNByRank[i-1];
      recvDataRankDisplacements[i] = recvDataRankDisplacements[i-1] + recvDataNByRank[i-1];
    }
  }
  /*
  if(rank == 0){
    std::cout << "Data per rank send, cum send, recv, cum recv\n";
    for(int i = 0; i < size; i++){
      std::cout << sendDataNByRank[i] << " " << sendDataRankDisplacements[i] << " " << recvDataNByRank[i] << " " << recvDataRankDisplacements[i] << "\n"; 
    }
  }
  */
  timerStart = std::clock();
  for(int i = 0; i < N; i++){
    microseconds = rand() % 20;
    usleep(microseconds);
    MPI_Alltoallv(&sendData[0],&sendDataNByRank[0],&sendDataRankDisplacements[0], MPI_INT,&recvData[0],&recvDataNByRank[0],&recvDataRankDisplacements[0], MPI_INT, MPI_COMM_WORLD);
  }
  if(rank==0){
    std::cout << "8 ranks, mid volume (" << totalVolume << ") each: " << (std::clock()-timerStart)/(double)(CLOCKS_PER_SEC) << std::endl;
  }

  for(int i = 0; i < size; i++){
    sendDataNByRank[i] = 0; 
    recvDataNByRank[i] = 0; 
    sendDataRankDisplacements[i] = 0;
    recvDataRankDisplacements[i] = 0;
    totalVolume += sendDataNByRank[i];
  }

  totalVolume = 0;
  sendDataNByRank[0] = 2; 
  recvDataNByRank[0] = 2; 
  for(int i = 0; i < size; i++){
    sendDataNByRank[i] = 2; 
    recvDataNByRank[i] = 2; 
    sendDataRankDisplacements[i] = sendDataRankDisplacements[i-1] + sendDataNByRank[i-1];
    recvDataRankDisplacements[i] = recvDataRankDisplacements[i-1] + recvDataNByRank[i-1];
    totalVolume += sendDataNByRank[i];
  }

  /*
  if(rank == 0){
    std::cout << "Data per rank send, cum send, recv, cum recv\n";
    for(int i = 0; i < size; i++){
      std::cout << sendDataNByRank[i] << " " << sendDataRankDisplacements[i] << " " << recvDataNByRank[i] << " " << recvDataRankDisplacements[i] << "\n"; 
    }
  }
  */

  timerStart = std::clock();
  for(int i = 0; i < N; i++){
    microseconds = rand() % 20;
    usleep(microseconds);
    MPI_Alltoallv(&sendData[0],&sendDataNByRank[0],&sendDataRankDisplacements[0], MPI_INT,&recvData[0],&recvDataNByRank[0],&recvDataRankDisplacements[0], MPI_INT, MPI_COMM_WORLD);
  }
  if(rank==0){
    std::cout << "All ranks, low volume (" << totalVolume << ") each: " << (std::clock()-timerStart)/(double)(CLOCKS_PER_SEC) << std::endl;
  }
 

  MPI_Finalize();

  std::cout << "Rank "<< rank << " is done!\n";
}
