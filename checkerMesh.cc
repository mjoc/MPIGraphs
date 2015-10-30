#include <iostream>
#include <vector>
#include <map>
#include <fstream>

int main(int argc, const char * argv[])
{
  //int nBigRows = 4;
  //int nBigCols = 6;
  //int nSmallRows = 5;
  //int nSmallCols = 5;
  //int nOverlapRows = 1;

  
  int nBigRows = 4;
  int nBigCols = 8;
  int nSmallRows = 10;
  int nSmallCols = 10;

  
  int iCurrentBlock = 0;

  
  std::ofstream adjFile, coordsFile, elFile; //, edgeListFile;
  //std::string dim = NumberToString(N);;

  char bigR[20];
  sprintf(bigR, "%d", nBigRows);
  
  char bigC[20];
  sprintf(bigC, "%d", nBigCols);
  
  char smallR[20];
  sprintf(smallR, "%d", nSmallRows);
  
  char smallC[20];
  sprintf(smallC, "%d", nSmallCols);
 
  std::string outDir("/home/users/mschpc/2013/oconnm28/Courses/GraphPartitioning/Code/Graphs/");
  std::string fileName(outDir + "checkerR" + bigR + "C" + bigC + "r" + smallR + "c"+ smallC + ".txt");
 
  adjFile.open(fileName.c_str());
  
  fileName= std::string(outDir + "checkerR" + bigR + "C" + bigC + "r" + smallR + "c" + smallC + "_coords.txt");
 
  coordsFile.open(fileName.c_str());

  fileName= std::string(outDir + "checkerR" + bigR + "C" + bigC + "r" + smallR + "c" + smallC + "_el.txt");
  elFile.open(fileName.c_str());

  std::vector<int> activeBlocks;
  int iBlock = 0;
  for(int iRow = 0 ; iRow < nBigRows; iRow++){
    iBlock = iRow * nBigCols - 1-iRow%2;
    for(int iCol = 0; iCol < nBigCols/2; iCol++){
      iBlock += 2;
      activeBlocks.push_back(iBlock);
      std::cout << iBlock << " ";
    }
  }
  std::cout << std::endl;

  
  //std::vector<int> activeBlocks= {1,3,5,6,8,10,13,15,17,18,20,22}; //{0,2,4,7,9,11,12,14,16,19,21,23};
  std::vector<bool> blockActive(nBigRows*nBigCols);
  for(std::vector<bool>::iterator it = blockActive.begin(); it != blockActive.end(); it++){
    *it = false;
  }
  for(std::vector<int>::iterator it = activeBlocks.begin(); it != activeBlocks.end(); it++){
    blockActive[*it] = true;
  }

  
  int nAdded = 0, iAddCol;
  std::map<std::pair<int,int>, int> vertices;
  std::map<int,std::pair<int,int> > vertices2;
  int iVertex = 0;
  for(int iRow = 0; iRow < (nBigRows*nSmallRows);iRow++){
    for(int iCol = 0; iCol < (nBigCols*nSmallCols);iCol++){
      iCurrentBlock = (iRow/nSmallRows * nBigCols) + (iCol/nSmallCols);
      if(blockActive[iCurrentBlock]){
        //std::cout << "Row " << iRow << "\n";

    
        if(((iCol+1) % nSmallCols) < 2 || ((iRow + 1) % nSmallRows < 2) ){
          if(((iCol+1) % nSmallCols) < 2 && ((iRow + 1) % nSmallRows < 2) ){
            if(iCol % nSmallCols == 0){
              if(iCol == 0){
                iAddCol = (nBigCols*nSmallCols)-1;
              }else{
                iAddCol = iCol-1;
              }
              std::cout << "Adding Row " << iRow << " Col " << iAddCol << "\n";
              nAdded++;
              vertices.insert(std::pair<std::pair<int,int>,int>(std::make_pair(iRow,iAddCol), iVertex));
              vertices2.insert(std::pair<int,std::pair<int,int> >(iVertex, std::make_pair(iRow,iAddCol)));
              iVertex++;
            }
            nAdded++;
            vertices.insert(std::pair<std::pair<int,int>,int>(std::make_pair(iRow,iCol), iVertex));
            vertices2.insert(std::pair<int,std::pair<int,int> >(iVertex, std::make_pair(iRow,iCol)));
            iVertex++;

            
            std::cout << "Row " << iRow << " Col " << iCol << " is corner point\n";
            if(iCol % nSmallCols == (nSmallCols-1)){
              iAddCol = (iCol+1)%(nBigCols*nSmallCols);
              std::cout << "Adding Row " << iRow << " Col " << iAddCol << "\n";
              nAdded++;
              vertices.insert(std::pair<std::pair<int,int>,int>(std::make_pair(iRow,iAddCol), iVertex));
              vertices2.insert(std::pair<int,std::pair<int,int> >(iVertex, std::make_pair(iRow,iAddCol)));
              iVertex++;
            }

            
          }else{
            std::cout << "Adding Row " << iRow << " Col " << iCol << "\n";
            nAdded++;
            vertices.insert(std::pair<std::pair<int,int>,int>(std::make_pair(iRow,iCol), iVertex));
            vertices2.insert(std::pair<int,std::pair<int,int> >(iVertex, std::make_pair(iRow,iCol)));
            iVertex++;
          }
        }else{
          std::cout << "Adding Row " << iRow << " Col " << iCol << "\n";
          nAdded++;
          vertices.insert(std::pair<std::pair<int,int>,int>(std::make_pair(iRow,iCol), iVertex));
          vertices2.insert(std::pair<int,std::pair<int,int> >(iVertex, std::make_pair(iRow,iCol)));
          iVertex++;
        }
      }
    }
    std::cout << "Row " << iRow << ", added so far: " << vertices.size() << std::endl;
    std::cout << "Row " << iRow << ", added so far: " << vertices2.size() << std::endl;
  }
  std::cout << "Added corners: " << nAdded << "\n";
  std::cout << "Added " <<  vertices.size() << " vertices, expected " << (nSmallRows* nSmallCols)*(nBigRows* nBigCols)/2 + (nBigRows* nBigCols)*2 << "\n";

  

  
  //for(std::map<std::pair<int,int>, int>::const_iterator it = vertices.begin(); it != vertices.end(); it++ ){
  //  std::cout << (*it).second << " (" << ((*it).first).first << ", " << ((*it).first).second << ")\n";
  //}

  

  
  //for(int iRow = 0; iRow < (nBigRows*nSmallRows);iRow++){
  //  for(int iCol = 0; iCol < (nBigCols*nSmallCols);iCol++){
  int iRow, iCol;
  iVertex = 0;
  for(std::map<int,std::pair<int,int> >::const_iterator it = vertices2.begin(); it != vertices2.end(); it++){
    iRow = ((*it).second).first;
    iCol = ((*it).second).second;
    if(vertices.find(std::make_pair(iRow,iCol)) != vertices.end()){

      
      int iCheck = 0;
      int iRowMinusOne = ((iRow + (nSmallRows* nBigRows))-1) % (nSmallRows* nBigRows);
      int iRowPlusOne = (iRow + 1)  % (nSmallRows* nBigRows);
      int iColMinusOne = (iCol + (nSmallCols* nBigCols) - 1) % (nSmallCols* nBigCols);
      int iColPlusOne = (iCol +1) %  (nSmallCols* nBigCols);

      
      //std::cout << "Checking (" << iRow << ", " << iCol << ")" << std::endl;
      adjFile << vertices.find(std::make_pair(iRow,iCol))->second;
      if(vertices.find(std::make_pair(iRowMinusOne,iCol)) != vertices.end()){
        //std::cout << "A ";
        adjFile << "; " << vertices.find(std::make_pair(iRowMinusOne,iCol))->second;
	elFile << vertices.find(std::make_pair(iRow,iCol))->second << " " << vertices.find(std::make_pair(iRowMinusOne,iCol))->second << "\n";
        iCheck++;
      }
      if(vertices.find(std::make_pair(iRow,iColPlusOne)) != vertices.end()){
        //std::cout << "B ";
        adjFile << "; " << vertices.find(std::make_pair(iRow,iColPlusOne))->second;
	elFile << vertices.find(std::make_pair(iRow,iCol))->second << " " << vertices.find(std::make_pair(iRow,iColPlusOne))->second << "\n";
        iCheck++;
      }
      if(vertices.find(std::make_pair(iRowPlusOne,iCol)) != vertices.end()){
        //std::cout << "C ";
        adjFile << "; " << vertices.find(std::make_pair(iRowPlusOne,iCol))->second;
	elFile << vertices.find(std::make_pair(iRow,iCol))->second << " " << vertices.find(std::make_pair(iRowPlusOne,iCol))->second << "\n";
        iCheck++;
      }
      if(vertices.find(std::make_pair(iRow,iColMinusOne)) != vertices.end()){
        //std::cout << "D ";
        adjFile << "; " << vertices.find(std::make_pair(iRow,iColMinusOne))->second;
	elFile << vertices.find(std::make_pair(iRow,iCol))->second << " " << vertices.find(std::make_pair(iRow,iColMinusOne))->second << "\n";
        iCheck++;
      }
      adjFile << std::endl;
      //std::cout << std::endl;
      //if(((iCol+1) % nSmallCols) > 1 && ((iRow + 1) % nSmallRows > 1) ){
      //  if(iCheck != 4){
      //    std::cout << " There are not 4 edge for vertex (" << iRow << ", " << iCol << ")" << std::endl;
      //  }
      //}
      iCurrentBlock = (iRow/nSmallRows * nBigCols) + (iCol/nSmallCols);
      if(!blockActive[iCurrentBlock]){
        std::cout << " There are " << iCheck << " edge for vertex (" << iRow << ", " << iCol << ")" << std::endl;
      }
    }else{
      std::cout << "@@@@@@@@@ERRRORRR@@@@@@@@@\n";
    }

 
    //std::cout << "check " <<  vertices.size() << " equals " << vertices2.size() << "\n";

 
    coordsFile << iVertex << " " << vertices2[iVertex].first << " " << vertices2[iVertex].second << std::endl;
    iVertex++;
  }

  

  
  adjFile.close();
  coordsFile.close();
  elFile.close();
  return 0;
}
