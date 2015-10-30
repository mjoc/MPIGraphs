#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>

//#define N 4


//#include <sstream>

template <typename T>
std::string NumberToString ( T Number )
{
	std::ostringstream ss;
	ss << Number;
	return ss.str();
}



int main(int argc, char **argv){
  
  if(argc != 2){
    std::cerr << "Please give integer argument\n";
    exit(1);
  }
  
  int N = atoi(argv[1]);
  
  std::ofstream adjFile, coordsFile;
  //std::string dim = NumberToString(N);;
  std::string fileName = std::string("torusMesh") + argv[1] + ".txt";
  std::string fileName2 = std::string("torusMesh") + argv[1] + "_coords.txt";
  adjFile.open(fileName.c_str());
  coordsFile.open(fileName2.c_str());
  int iVertex = 0;
  for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++){
      int k = (i*N)+j ;
      adjFile << k;
      
      if((k - N)>-1){
        adjFile << ";"<< k - N;
      }else{
        adjFile << ";"<< N*N + (k - N);
      }
      
      if((k + 1) < ((k/N)+1)*N){
        adjFile << ";"<< k + 1 ;
      }else{
        adjFile << ";" << ((k/N)*N) + (k+1)%N;
      }
      
      if((k + N) < N*N){
        adjFile << ";"<< k + N ;
      }else{
        adjFile << ";"<< (k + N)%N;
      }
      
      if((k-1) >= (k/N)*N ){
        adjFile << ";"<< k - 1 ;
      }else{
        adjFile << ";"<<  ((k/N)+1)*N -1;
      }
      
      adjFile << std::endl;
      coordsFile << iVertex << " " << i << " " << j << std::endl;
      iVertex++;
    }
  }

  adjFile.close();
  coordsFile.close();
  std::cout << "Data written to " << fileName << std::endl;
	return 0;
}
