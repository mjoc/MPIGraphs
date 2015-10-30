#include "GraphPartitioning.hpp"
#include "Sparse.hpp"
#include "RngWrapper.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <cmath>
#include <ctime>
#include <numeric>
#include <algorithm>


#define EPSILON 0.000001

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}


void GraphPartitioning::normaliseL2(std::vector<double>& vect){
  double norm = normL2(vect);
  if(norm> EPSILON){
    for(std::vector<double>::iterator it = vect.begin(); it < vect.end(); it++){
      *it /= norm;
    }
  }
}

double GraphPartitioning::normL2(const std::vector<double> vect){
  double norm = 0.0;
  for(std::vector<double>::const_iterator it = vect.begin(); it < vect.end(); it++){
    norm += (*it) * (*it);
  }
  return sqrt(norm);
}

std::vector<double> GraphPartitioning::lanczos(const Sparse<double>& A, std::vector<double>& diag, std::vector<double>& offDiag, RngWrapper& random, std::vector<double> rvs, Matrix<double>& Q){
  
  std::clock_t  mainStart, localStart;
  mainStart = std::clock();
  localStart = std::clock();
  size_t dim = A.nRows();
  
  std::vector<std::vector<double> > Qkeep(A._nCols);
  
  bool randomSupplied = false;
  if(rvs.size()!=0){
    if(rvs.size() != dim){
      throw std::logic_error("Error with supplied random variables");
    }else{
      randomSupplied= true;
    }
  }
  std::vector<double> qk(dim);

  diag.resize(0);
  offDiag.resize(0);
  
  std::vector<double> rk(dim);
  if(randomSupplied){
    rk = rvs; //random.getUniformRandom();
  }else{
    rk = random.getUniform01Randoms(dim);
  }
  std::vector<double> temp(rk);
  std::vector<double> orthoCorrection(dim);

  normaliseL2(rk);
  std::vector<double> Aqk;
  
  for(size_t k = 0; k < dim; k++){
    

    if(k %100 == 0){
      
      std::cout << "Lanczos iteration " << k << " of " << dim << ", previous took "  << (std::clock() - localStart) / (double)(CLOCKS_PER_SEC) << " seconds" <<std::endl;
      localStart = std::clock();
    }
    if( k == 0){
      qk = rk;
    }else{
      qk = scalerDiv(rk,offDiag[k-1]);
    }
    Qkeep[k] = qk;
    Aqk = mult(A,qk);
    if(k == 0){
      rk = Aqk;
    }else{
      rk = diff(Aqk,scalerMult(Qkeep[k-1],offDiag[k-1]));
    }
    
    diag.push_back(dot(qk,rk));
    rk = diff(rk,scalerMult(qk,diag[k]));
    
    for(size_t i = 0; i <= k; i++){
      if(i ==0){
        orthoCorrection = scalerMult(Qkeep[i], dot(rk, Qkeep[i]));
      }else{
        std::vector<double> check = scalerMult(Qkeep[i], dot(rk, Qkeep[i]));
        orthoCorrection =  sum(orthoCorrection, check);
      }
    }
    rk = diff(rk, orthoCorrection);
    
    
    for(size_t i = 0; i <= k; i++){
      if(i ==0){
        orthoCorrection = scalerMult(Qkeep[i], dot(rk, Qkeep[i]));
      }else{
        std::vector<double> check = scalerMult(Qkeep[i], dot(rk, Qkeep[i]));
        orthoCorrection =  sum(orthoCorrection, check);
      }
    }
    rk = diff(rk, orthoCorrection);
    
    offDiag.push_back(normL2(rk));
    //Qkeep[k] = qk;
    
    //k++;
    if(k == dim-1){
      //keepGoing = false;
      offDiag.pop_back();
    }
  }
  //std::cout << "Lanczos, fixing up Q\n";
  Q.resize(A._nRows,A._nCols);
  for(int i= 0;i < A._nRows; i++){
    for(int j = 0; j < A._nCols; j++){
      Q._data[i*A._nCols+j] = Qkeep[j][i];
    }
  }
  //std::cout << "Lanczos, Q fixed up\n";
  std::cout << "Lanczos finished in " << (std::clock() - mainStart) / (double)(CLOCKS_PER_SEC) << " seconds" << std::endl;
  return temp;
}




void GraphPartitioning::getEigVals(const std::vector<double> &diag, const std::vector<double> &offDiag,const int n , std::vector<double>& eigVals){
  eigVals.resize(n);
  std::vector<double> brackets;
  eigValsBrackets(diag, offDiag, n, brackets);
  for(int i = 0; i < n; i ++){
    std::cout << "Looking for eigen values between (" << brackets[i] << " , " << brackets[i+1] << ")\n";
    try{
      eigVals[i] = bisection(diag,offDiag,brackets[i],brackets[i+1]);
      std::cout << "eigVals[" << i <<"] = " <<  eigVals[i]  << std::endl;
    }catch(std::logic_error& err){
      std::cerr << err.what();
    }
  }
  return;
}

double GraphPartitioning::getSecondEigenValue(const std::vector<double> &diag, const std::vector<double> &offDiag, double& result){
  bool allOk = true;
  double lowerBracket, upperBracket;

  
  std::vector<double>brackets(2);
  eigValsBrackets(diag, offDiag, 2, brackets);
  //secondEigValueBrackets(diag, offDiag, lowerBracket,upperBracket);
  lowerBracket = brackets[1];
  upperBracket = brackets[2];
  std::cout << "Looking for eigen value between (" << lowerBracket << " , " << upperBracket << ")\n";

  result = bisection(diag,offDiag,lowerBracket,upperBracket);
  std::cout << "eigValue: " << result  << std::endl;
  return allOk;
}

void GraphPartitioning::eigValsBrackets(const std::vector<double> &diag, const std::vector<double>& offDiag,const int n, std::vector<double>& brackets){
  
  double lb = 0.0, ub = 0.0;
  double eigVal, step;
  size_t maxRecursion = 100;
  size_t nEigVals;
  eigValBounds(diag,offDiag,lb,ub);
  brackets.resize(n+1);
  brackets[0] = lb;
  
  for(size_t i = n; i > 0; i--){
    eigVal = (lb+ub)/2.0;
    step = (ub-lb)/2.0;
    for(size_t j = 0; j <  maxRecursion; j++){
      nEigVals = countEigVals(diag, offDiag, eigVal);
      step = step/2;
      if(nEigVals < i){
        eigVal = eigVal + step;
      }else if(nEigVals > i){
        eigVal = eigVal - step;
      }else{
        break;
      }
    }
    brackets[i] = eigVal;
  }
  
}

int GraphPartitioning::countEigVals(const std::vector<double> &diag, const std::vector<double> &offDiag, const double lambda){
  int count, sign,oldSign;
  size_t n;
  std::vector<double>sturms;
  
  sturmSeq(diag, offDiag, lambda, sturms);
  oldSign = (int)sgn(sturms[0]);
  n = sturms.size();
  count = 0;
  for(size_t i = 1; i < n; i++){
    if(sturms[i] != sturms[i]){
      //Break on NaN
      break;
    }
    sign = sgn(sturms[i]);
    if(fabs(sign) < EPSILON) sign = -1*oldSign;
    if( sign*oldSign < 0) count++;
    oldSign = sign;
  }
  return count;
}


// Gerschgorin
void GraphPartitioning::eigValBounds(const std::vector<double> &diag, const std::vector<double> &offDiag, double &lower, double &upper){
  size_t n = diag.size();
  
  lower = diag[0] - fabs(offDiag[0]);
  upper = diag[0] + fabs(offDiag[0]);
  double eig = 0;
  for(size_t i = 1; i < n-1; i++){
    eig = diag[i] - fabs(offDiag[i]) - fabs(offDiag[i-1]);
    lower = eig < lower ? eig : lower;
    eig = diag[i] + fabs(offDiag[i]) + fabs(offDiag[i-1]);
    upper = eig > upper ? eig : upper;
  }
  eig = diag[n-1] - fabs(offDiag[n-2]);
  lower = eig < lower ? eig : lower;
  eig = diag[n-1] + fabs(offDiag[n-2]);
  upper = eig > upper ? eig : upper;
  return;
}

void GraphPartitioning::sturmSeq(const std::vector<double> &diag, const std::vector<double> &offDiag, const double lambda, std::vector<double> &seq){
  size_t n = diag.size()+1;
  seq.resize(n);
  seq[0] = 1;
  seq[1] = diag[0] - lambda;
  for(size_t i = 2; i < n; i++){
    seq[i] = (diag[i-1] - lambda)*seq[i-1] - (offDiag[i-2]*offDiag[i-2])*seq[i-2];
  }
  return;
}

double GraphPartitioning::lastSturm(const std::vector<double> &diag, const std::vector<double> &offDiag, double lambda){
  std::vector<double> sturms;
  sturmSeq(diag,offDiag,lambda,sturms);
  
  return sturms[sturms.size()-1];
}


double GraphPartitioning::bisection(const std::vector<double> &diag, const std::vector<double> &offDiag, double lower, double upper){
  double fLower = lastSturm(diag,offDiag,lower);
  double fUpper = lastSturm(diag,offDiag,upper);
  double estimatedRoot;
  
  //std::cout << "Bisection: looking between " << lower << " and " << upper << std::endl;
  //std::cout << "fLower " << fLower << " and fUpper " << fUpper << std::endl;
  
  if(fabs(upper-lower) < EPSILON){
    estimatedRoot = lower + 0.5*(upper-lower);
  }else{
    if(fabs(fLower) < EPSILON){
      estimatedRoot = lower;
    }else if(fabs(fUpper) < EPSILON){
      estimatedRoot = upper;
    }else{
      if(fLower * fUpper > 0.0){
        std::cout << "The supplied bounds evaluate to the same sign, can't do bisection\n";
        estimatedRoot = lower + (upper-lower)/2;
      }
    
      double split = lower + 0.5*(upper-lower);
      double fSplit = lastSturm(diag,offDiag,split);
      if(fLower * fSplit < 0.0){
        upper = split;
      }else{
        lower = split;
      }
    estimatedRoot = bisection(diag, offDiag, lower, upper);
    }
  }
  return estimatedRoot;
}



bool GraphPartitioning::inversePower(Sparse<double>& mat, double eigenEstimate, std::vector<double>& eigenVector, double& eigenValue, RngWrapper& rng){
  
  
  double epsilon = 0.00000001;
  
  bool converged = true;
  std::vector<int> dataEntryAdded;
  int iData = 0;
  std::vector<int> diagonals;
  
  eigenVector = rng.getUniform01Randoms(mat._nRows);
  
  std::map<std::pair<int,int>, double> matrix;
  bool diagFound = false;
  for(int iRow = 0; iRow < mat._nRows;iRow++){
    iData = mat._rowIndex[iRow];
    diagFound = false;
    while(iData < mat._rowIndex[iRow+1]){
      if(mat._cols[iData] == iRow){
        matrix[std::make_pair(iRow,iRow)] = mat._data[iData]-eigenEstimate;
        diagFound = true;
      }else{
        matrix[std::make_pair(iRow,mat._cols[iData])] = mat._data[iData];
      }
      iData++;
    }
    if(!diagFound){
      matrix[std::make_pair(iRow,iRow)] = -eigenEstimate;
    }
  }
  
  std::vector<int> rowOrder(mat._nRows);
  for(int i = 0; i < (int)rowOrder.size();i++){
    rowOrder[i] = i;
  }
  bool allOk = true;
  allOk = luDecomp(matrix, mat._nRows,  mat._nRows,  rowOrder);
  if(!allOk){
    std::cout << "Problem with luDecomp\n";
  }
  
  bool notConverged = true;
  int nIter = 0, maxIter = 100;
  double maxAbsX, lambda = 0, dc, dv;
  size_t iMaxAbsX;
  double bigX;
  std::vector<double> eigVecCalc(mat._nRows);
  while (nIter < maxIter && notConverged && allOk){
    allOk = luSolve(matrix, mat._nRows, mat._nRows, rowOrder, eigenVector, eigVecCalc);
    if(allOk){
      maxAbsX = fabs(eigVecCalc[0]);
      iMaxAbsX = 0;
      for(size_t i = 1; i < eigVecCalc.size(); i++){
        if(eigVecCalc[i]>maxAbsX){
          maxAbsX = fabs(eigVecCalc[i]);
          iMaxAbsX = i;
        }
      }
      bigX = eigVecCalc[iMaxAbsX];
      dc = fabs(lambda - bigX);
      scaleVector(eigVecCalc,1/bigX);
      dv = normL2(diff(eigenVector,eigVecCalc));
      eigenVector = eigVecCalc;
      lambda = bigX;
      if(fmax(dc,dv) < epsilon){
        notConverged = false;
      }
      lambda = eigenEstimate + 1/bigX;
      nIter++;
    }else{
      std::cout << " Problem with LUsolve\n";
    }
  }
  if(converged){
    eigenValue = lambda;
  }
  return converged;
}

bool GraphPartitioning::luDecomp(std::map<std::pair<int,int>, double>& matrix, int nRows, int nCols, std::vector<int>& rowOrder){
  bool allOk = true;
  double max = 0.0;
  int iPivot;
  double det1 = 1.0, dataValue;
   
  for(int iRow = 0; iRow < nRows-1;iRow++){
    
    if(matrix.find(std::make_pair(iRow,iRow)) == matrix.end()){
      max = 0.0;
      iPivot = iRow;
    }else{
      max = fabs(matrix.find(std::make_pair(iRow,iRow))->second);
      iPivot = iRow;
    }
   
    for(int jRow = iRow+1; jRow < nRows; jRow++){
      if(fabs(matrix.find(std::make_pair(jRow,iRow))->second) > max){
        max = fabs(matrix.find(std::make_pair(jRow,iRow))->second);
        iPivot = jRow;
      }
    }
    if(iPivot > iRow){
      int temp = rowOrder[iRow];
      rowOrder[iRow] = rowOrder[iPivot];
      rowOrder[ iPivot] = temp;
      det1 *= -1.0;
    }
    if(matrix.find(std::make_pair(rowOrder[iRow],iRow)) == matrix.end()){
      dataValue = 0.0;
    }else{
      dataValue = matrix.find(std::make_pair(rowOrder[iRow],iRow))->second;
    }
    det1 = det1*dataValue;
    if(fabs(det1) < 0.00000001){
      allOk = false;
      std::cout << "determinant is zero, stopping\n";
      break;
    }
    double factor;
    for(int jRow = iRow+1; jRow < nRows; jRow++){
      
      if(matrix.find(std::make_pair(rowOrder[jRow],iRow)) == matrix.end()){
        dataValue = 0.0;
      }else{
        dataValue = matrix.find(std::make_pair(rowOrder[jRow],iRow))->second;
      }
      if(matrix.find(std::make_pair(rowOrder[iRow],iRow)) == matrix.end()){
        allOk = false;
        std::cout << "Denominator: ("<< rowOrder[iRow] << ", " << iRow << ") is zero, stopping\n";
        break;
      }else{
        factor = dataValue/matrix.find(std::make_pair(rowOrder[iRow],iRow))->second;
      }
  
      matrix[std::make_pair(rowOrder[jRow],iRow)] =  factor;
       
      for(int iCol = iRow + 1; iCol < nRows; iCol++){
        if(matrix.find(std::make_pair(rowOrder[jRow],iCol)) == matrix.end()){

          matrix[std::make_pair(rowOrder[jRow],iCol)] = 0.0;
        }
        
        if(matrix.find(std::make_pair(rowOrder[iRow],iCol)) == matrix.end()){
          dataValue = 0.0;
        }else{
          dataValue = matrix.find(std::make_pair(rowOrder[iRow],iCol))->second;
        }
        
        matrix[std::make_pair(rowOrder[jRow],iCol)] -= factor * dataValue;
      }
    }
  }
  
  return allOk;
}

bool GraphPartitioning::luSolve(std::map<std::pair<int,int>, double> A, int nRows, int nCols, std::vector<int> rowOrder, std::vector<double>& b, std::vector<double>& X){
  bool allOk = true;

  
  
  std::vector<double> Y(nRows);
  double dataValue = 0.0 , sumProduct = 0.0;
  Y[0] = b[rowOrder[0]];
  for(int i = 1; i < nRows;i++){
    sumProduct = 0.0;
    for(int j = 0; j < i; j++){
      if(A.find(std::make_pair(rowOrder[i],j)) == A.end()){
        std::cout << "Cannot find (" << rowOrder[i] << ", " << j << ")\n";
        dataValue = 0.0;
      }else{
        dataValue = A.find(std::make_pair(rowOrder[i],j))->second;
      }
      sumProduct += dataValue * Y[j];
    }
    
    Y[i] = b[rowOrder[i]] - sumProduct;
  }
  
  if(A.find(std::make_pair(rowOrder[nRows-1],nRows-1)) == A.end()){
    allOk =false;
  }else{
    dataValue =  A.find(std::make_pair(rowOrder[nRows-1],nRows-1))->second;
  }

  if(allOk){
    X[nRows-1] = Y[nRows-1]/dataValue;
    for(int i = nRows - 2; i > -1; i--){
      X[i] = Y[i];
      sumProduct = 0.0;
      for(int j = i+1; j < nRows; j++){
        if(A.find(std::make_pair(rowOrder[i],j)) == A.end()){
          dataValue = 0.0;
        }else{
          dataValue =  A.find(std::make_pair(rowOrder[i],j))->second;
        }
        X[i] -= dataValue * X[j];
      }
      
      if(A.find(std::make_pair(rowOrder[i],i)) == A.end()){
        allOk = false;
        break;
      }else{
        dataValue =  A.find(std::make_pair(rowOrder[i],i))->second;
      }
      
      if(allOk){
        X[i] /= dataValue;
      }
      
    }
    
    
  }

  return allOk;
}



bool GraphPartitioning::inversePowerTri(std::vector<double> diag, std::vector<double>& lowerDiag, double eigenEst, std::vector<double>& eigenVector, double& eigenValue, RngWrapper& rng){
  bool converged = false;
  double tol = 0.00000001;
  int matIter = 100;
  
  
  std::vector<double> upperDiag = lowerDiag;
  scalerSubtract(diag,eigenEst);

  luDecompTri(diag,lowerDiag,upperDiag);
  

  eigenVector = rng.getUniform01Randoms(diag.size());
  normaliseL2(eigenVector);
  std::vector<double> oldEig;
  int xSign = 0;
  double eigNorm;
   
  for(int i = 0; i < matIter; i++){
 
    oldEig = eigenVector;
    
 
    luSolveTri(diag, lowerDiag, upperDiag, eigenVector);
      
    eigNorm = normL2(eigenVector);
    normaliseL2(eigenVector);
      
    xSign = dot(oldEig,eigenVector) < 0?-1:1;
    if(xSign==-1){
      scaleVector(eigenVector,-1.0);
    }
      
    if(normL2(diff(oldEig,eigenVector)) < tol){
      std::cout << "Inverse Power Converged, norm:" << eigNorm << "\n";
      eigenValue = eigenEst + (double)xSign/eigNorm;
      converged = true;
      break;
    }
    
  }
  return converged;
}

void GraphPartitioning::luDecompTri(std::vector<double>& diag, std::vector<double>& lowerDiag,const std::vector<double>& upperDiag){
  
  for(size_t i = 1; i < diag.size();i++){
    lowerDiag[i-1] /= diag[i-1];
    diag[i] -= lowerDiag[i-1] * upperDiag[i-1];
  }
  return;
}

void GraphPartitioning::luSolveTri(const std::vector<double>& diag, const std::vector<double>& lowerDiag, const std::vector<double>& upperDiag, std::vector<double>& b){
  
  for(size_t i = 1; i < b.size();i++){
    b[i] = b[i] - lowerDiag[i-1]*b[i-1];
  }
  
  b[diag.size()-1] /= diag[diag.size()-1];
  
  for(int i = (int)diag.size()-2; i > -1;i--){
    b[i] = (b[i]-upperDiag[i]*b[i+1])/diag[i];
  }
}





double GraphPartitioning::dot(const std::vector<double>& a, const std::vector<double>& b){
  if(a.size() != b.size()){
    throw std::logic_error("Vectors of unequal length in dot!!");
  }

  return std::inner_product(a.begin(),a.end(),b.begin(),0.0);
}


std::vector<double> GraphPartitioning::sum(const std::vector<double>& a, const std::vector<double>& b){
  if(a.size() != b.size()){
    throw std::logic_error("Vectors of unequal length in dot!!");
  }

  std::vector<double> result;
  result.reserve(a.size());
  
  std::transform(a.begin(), a.end(), b.begin(),
                 std::back_inserter(result), std::plus<double>());
  return result;
}

std::vector<double> GraphPartitioning::diff(const std::vector<double>& a, const std::vector<double>& b){
  if(a.size() != b.size()){
    throw std::logic_error("Vectors of unequal length in dot!!");
  }
  
  std::vector<double> result;
  result.reserve(a.size());
  
  std::transform(a.begin(), a.end(), b.begin(),
                 std::back_inserter(result), std::minus<double>());
  return result;
}




std::vector<double> GraphPartitioning::scalerDiv(const std::vector<double>& vec, const double denom){
  std::vector<double> result(vec);
  for(std::vector<double>::iterator it = result.begin(); it < result.end(); ++it){
    *it /= denom;
  }
  return result;
  
}
std::vector<double> GraphPartitioning::mult(const Sparse<double>& mat,const std::vector<double>& vect){
  
  if(mat._nCols != (int)vect.size()){
    throw std::logic_error("Vectors of unequal length in dot!!");
  }
  
  
  std::vector<double> result(mat._nRows);
  double sumproduct;
  int iDataIndex, iCol;
  for(int iRow = 0; iRow < mat._nRows; iRow++){
    if(mat._rowIndex[iRow] < mat._rowIndex[iRow+1]){
      sumproduct = 0.0;
      for(iDataIndex = mat._rowIndex[iRow]; iDataIndex< mat._rowIndex[iRow+1];  iDataIndex++){
        iCol = mat._cols[iDataIndex];
        sumproduct += (double)mat._data[iDataIndex] * vect[iCol];

      }
      result[iRow] = sumproduct;
    }
  }
  return result;
}

std::vector<double> GraphPartitioning::mult(const Matrix<double>& mat,const std::vector<double>& vect){
  
  if(mat._nCols != (int)vect.size()){
    throw std::logic_error("Vectors of unequal length in dot!!");
  }
  
  std::vector<double> result(mat._nRows);
  for(int iRow = 0; iRow < mat._nRows; iRow++){
    result[iRow] = 0.0;
    for(int iCol = 0; iCol < mat._nCols;  iCol++){
      result[iRow] += mat._data[iRow*mat._nCols + iCol]* vect[iCol];
    }
  }
  return result;
}

std::vector<double> GraphPartitioning::scalerMult(const std::vector<double>& vec, double scaler){
  std::vector<double> result(vec);
  for(std::vector<double>::iterator it = result.begin(); it < result.end(); ++it){
    *it *= scaler;
  }
  return result;
}

void GraphPartitioning::scaleVector(std::vector<double>& vec, double scaler){
  for(std::vector<double>::iterator it = vec.begin(); it < vec.end(); ++it){
    *it *= scaler;
  }
}

void GraphPartitioning::scalerSubtract(std::vector<double>& vec, double scaler){
  for(std::vector<double>::iterator it = vec.begin(); it < vec.end(); ++it){
    *it -= scaler;
  }
  return;
}


bool GraphPartitioning::partitionGraph(const Sparse<double>& laplacian, std::vector<int>& groups, int nSplits, int useMedian, RngWrapper& rng){
  bool allOk = true;
  char message[200];
  std::vector<double> offDiag;
  std::vector<double> diag;
  std::vector<double> dummy;
  Matrix<double> Qmat;
  std::clock_t  mainStart, localStart;
  
  mainStart = std::clock();
  
  sprintf (message,  "Lanczos calculation");
  _runlog->writelog(message);
  
  lanczos(laplacian,diag,offDiag, rng, dummy, Qmat);
  std::cout << "Lanczos done\n";
  sprintf (message,  "Time of calculation Lanczos of size %d: %f seconds", laplacian._nRows, (std::clock() - mainStart) / (double)(CLOCKS_PER_SEC));
  _runlog->writelog(message);
  
  localStart = std::clock();
  
  double lb, ub;
  eigValBounds(diag,offDiag,lb,ub);
  std::cout << "Eig Val bounds done\n";
  sprintf (message,  "Time of calculation eigValBounds of size %d: %f seconds", laplacian._nRows, (std::clock() - localStart) / (double)(CLOCKS_PER_SEC));
  _runlog->writelog(message);
  
  localStart = std::clock();
  double lanczosEigenvalue;
  sprintf (message,  "Looking for second eigenvalue");
  _runlog->writelog(message);
  
  allOk = getSecondEigenValue(diag, offDiag, lanczosEigenvalue);
  std::cout << "Second Eigen bounds done\n";
  if(!allOk){
    std::cerr << "Calculating second eigenvalue by Lanczos failed!\n";
    return 1;
  }
  
  sprintf (message,  "Time of calculation Second Eigenvalue of size %d: %f seconds", laplacian._nRows, (std::clock() - localStart) / (double)(CLOCKS_PER_SEC));
  _runlog->writelog(message);
  
  bool gotEigenVec = false;
  double eigenValueRecalc;
  std::vector<double> eigenVectorTri(diag.size());
  localStart = std::clock();
  gotEigenVec = inversePowerTri(diag, offDiag, lanczosEigenvalue, eigenVectorTri , eigenValueRecalc, rng);
  sprintf (message,  "Time of calculation inversePowerTri of size %d: %f seconds", laplacian._nRows, (std::clock() - localStart) / (double)(CLOCKS_PER_SEC));
  _runlog->writelog(message);
  
  if(gotEigenVec){
    std::cout<< "tridiagonal Inverse Power converged to T eigenvector\n";
  }else{
    std::cout << "Something wrong with tridiagonal Inverse Power\n";
  }
  
  std::cout << "Multiplying Q * T-eigenvector to get Laplacian Eigenvector \n";
  localStart = std::clock();
  std::vector<double> ritz = mult(Qmat,eigenVectorTri);
  sprintf (message,  "Time of calculation Q*T-eigenvector of size %d: %f seconds", laplacian._nRows, (std::clock() - localStart) / (double)(CLOCKS_PER_SEC));
  _runlog->writelog(message);
  
  localStart = std::clock();
  groups.resize(ritz.size());
  
  std::vector<int> rowIndexA;
  std::vector<double> dataA;
  std::vector<int> colsA;
  std::vector<int> verticesA;

  double cutoff = 0.0;
  if(useMedian==1){
    sprintf (message,  "Using median of Fiedler Vector to split");
    _runlog->writelog(message);

    std::vector<double> ritz2(ritz);
    sort(ritz2.begin(), ritz2.end());

    
    int middle = (int)ritz2.size()/2;
    cutoff = ritz2.size() % 2 == 0 ? (ritz2[middle] + ritz2[middle-1]) / 2 : ritz2[middle];
  }
  
  for(int i = 0; i < (int)ritz.size(); i++){
    if(ritz[i] < cutoff){
      groups[i] = 0;
    }else{
      groups[i] = 1;
    }
  }
  
  if(nSplits > 1){
    for(int iSplit = 1; iSplit < nSplits; iSplit++){
      for(int iPartition = 0; iPartition < (int)pow(2,iSplit);iPartition++){
        std::cout << "Subdivision: " << iSplit << " partition " << iPartition << " of " << (int)pow(2,iSplit) << "\n";
        verticesA.resize(0);
        
        int nGroupSize = 0;
        for(int iVertex = 0;iVertex < (int)groups.size();iVertex++){
          if(groups[iVertex] == iPartition){
            verticesA.push_back(iVertex);
            nGroupSize++;
          }
        }

        rowIndexA.resize(verticesA.size()+1);
        dataA.resize(0);
        colsA.resize(0);

        std::vector<int> rowCols;
        rowIndexA[0] = 0;
        int iVertexA = 0;
        localStart = std::clock();
        std::vector<int>::const_iterator findIterator;
        for(std::vector<int>::iterator it = verticesA.begin();it  != verticesA.end();++it){
          rowCols.resize(0);
          for(int iVertex2 = laplacian._rowIndex[*it] ; iVertex2 < laplacian._rowIndex[(*it)+1]; iVertex2++){
            if(*it != laplacian._cols[iVertex2]){
              if(groups[laplacian._cols[iVertex2]] ==iPartition){
                findIterator = find(verticesA.begin(),verticesA.end(),laplacian._cols[iVertex2]);
                rowCols.push_back((int)(findIterator - verticesA.begin()));
              }
            }
          }
          bool diagEntered = false;
          for(int iVertex2 = 0; iVertex2 < (int)rowCols.size(); iVertex2++){
            if(!diagEntered && rowCols[iVertex2]> iVertexA){
              dataA.push_back(rowCols.size());
              colsA.push_back(iVertexA);
              diagEntered = true;
            }
            dataA.push_back(-1);
            colsA.push_back(rowCols[iVertex2]);
          }
          if(!diagEntered){
            dataA.push_back(rowCols.size());
            colsA.push_back(iVertexA);
          }
          rowIndexA[iVertexA+1] = (int)rowIndexA[iVertexA] + (int)rowCols.size() + 1;
          iVertexA++;
        }
        
        sprintf (message,  "Time of reorganising data size %d: %f seconds", laplacian._nRows, (std::clock() - localStart) / (double)(CLOCKS_PER_SEC));
        _runlog->writelog(message);
        {
          Sparse<double> laplacianA;
          double lb, ub;
        
          std::vector<int> groupsA(verticesA.size());
        
          if(dataA.size()>0){
	    std::cout << "Size " << verticesA.size() << std::endl;
            laplacianA = Sparse<double>((int)verticesA.size(),(int)verticesA.size(),dataA,colsA,rowIndexA);
          
            sprintf (message,  "Lanczos calculation");
            _runlog->writelog(message);
          
            lanczos(laplacianA,diag,offDiag, rng, dummy, Qmat);
            std::cout << "Lanczos done\n";
            sprintf (message,  "Time of calculation Lanczos of size %d: %f seconds", laplacianA._nRows, (std::clock() - mainStart) / (double)(CLOCKS_PER_SEC));
            _runlog->writelog(message);
          
            localStart = std::clock();
          
            eigValBounds(diag,offDiag,lb,ub);
            std::cout << "Eig Val bounds done\n";
            sprintf (message,  "Time of calculation eigValBounds of size %d: %f seconds", laplacianA._nRows, (std::clock() - localStart) / (double)(CLOCKS_PER_SEC));
            _runlog->writelog(message);
          
            localStart = std::clock();
            double lanczosEigenvalue;
            sprintf (message,  "Looking for second eigenvalue");
            _runlog->writelog(message);
          
            allOk = getSecondEigenValue(diag, offDiag, lanczosEigenvalue);
            std::cout << "Second Eigen bounds done\n";
            if(!allOk){
              std::cerr << "Calculating second eigenvalue by Lanczos failed!\n";
              return 1;
            }
          
            sprintf (message,  "Time of calculation Second Eigenvalue of size %d: %f seconds", laplacian._nRows, (std::clock() - localStart) / (double)(CLOCKS_PER_SEC));
            _runlog->writelog(message);
          
            bool gotEigenVec = false;
            double eigenValueRecalc;
            std::vector<double> eigenVectorTri(diag.size());
            localStart = std::clock();
            gotEigenVec = inversePowerTri(diag, offDiag, lanczosEigenvalue, eigenVectorTri , eigenValueRecalc, rng);
            sprintf (message,  "Time of calculation inversePowerTri of size %d: %f seconds", laplacian._nRows, (std::clock() - localStart) / (double)(CLOCKS_PER_SEC));
            _runlog->writelog(message);
          
            if(gotEigenVec){
              std::cout<< "tridiagonal Inverse Power converged to T eigenvector\n";
            }else{
              std::cout << "Something wrong with tridiagonal Inverse Power\n";
            }
          
            std::cout << "Multiplying Q * T-eigenvector to get Laplacian Eigenvector \n";
            localStart = std::clock();
            std::vector<double> ritz = mult(Qmat,eigenVectorTri);
            sprintf (message,  "Time of calculation Q*T-eigenvector of size %d: %f seconds", laplacian._nRows, (std::clock() - localStart) / (double)(CLOCKS_PER_SEC));
            _runlog->writelog(message);
	    
          
            if(allOk){
	      std::vector<double> ritz2(ritz);
              sort(ritz2.begin(), ritz2.end());
	      double cutoff = 0.0;
              if(useMedian==1){
		int middle = (int)ritz2.size()/2;
		cutoff = ritz2.size() % 2 == 0 ? (ritz2[middle] + ritz2[middle-1]) / 2 : ritz2[middle];
	      }
	      for(int i = 0; i < (int)ritz.size();i++){
		if(ritz[i] > cutoff){
		  groups[verticesA[i]] += (int)pow(2,iSplit);
		}
	      }
            }
          }
        }
      }
    }
  }
  sprintf (message,  "Time of full split graph function at size %d: %f seconds",  laplacian._nRows,(std::clock() - mainStart) / (double)(CLOCKS_PER_SEC));
  _runlog->writelog(message);

  
  return allOk;

  }
