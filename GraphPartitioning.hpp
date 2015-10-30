#ifndef __MOC_Eigen_
#define __MOC_Eigen_
#include <vector>
#include <map>

//This class includes a bunch of functions created along the way
//which are not used in the final project, maining full matrix algorithms 
//where tri diag functions were used instead

class Runlog;

template<typename DataType>
class Sparse;

class RngWrapper;

template<typename DataType>
struct Matrix
{
  // should probably be hidden away, and the class would
  // provide `at` and `operator()` for access
  int _nCols, _nRows;
  std::vector<DataType> _data;
  
  Matrix(){
    _nRows = 0;
    _nCols = 0;
  }
  Matrix(int nRows, int nCols){
    _nRows = nRows;
    _nCols = nCols;
    _data.resize(nRows * nCols);
  }
  void resize(int nRows, int nCols){
    _nRows = nRows;
    _nCols = nCols;
    _data.resize(nRows * nCols);
  }
  
};


class GraphPartitioning{
  Runlog *_runlog;
public:
  GraphPartitioning(Runlog *runlog){_runlog = runlog;};
  
  std::vector<double> lanczos(const Sparse<double>& , std::vector<double>& , std::vector<double>& , RngWrapper&, std::vector<double>, Matrix<double>& Q2);

  void eigValsBrackets(const std::vector<double> &alpha, const std::vector<double> &beta, const int n, std::vector<double>& brackets);

  void getEigVals(const std::vector<double> &alpha, const std::vector<double> &beta,const int, std::vector<double>& eigVals);
  
  double getSecondEigenValue(const std::vector<double> &diag, const std::vector<double> &offDiag, double& result);

  void eigValBounds(const std::vector<double> &alpha, const std::vector<double> &beta, double &lower, double &upper);
  
  void secondEigValueBrackets(const std::vector<double> &diag, const std::vector<double>& offDiag, double& lowerBrackets, double& upperBracket );

  void sturmSeq(const std::vector<double> &alpha, const std::vector<double> &beta, const double lambda, std::vector<double> &seq);

  double lastSturm(const std::vector<double> &alpha, const std::vector<double> &beta, double lambda);

  int countEigVals(const std::vector<double> &alpha, const std::vector<double> &beta, const double lambda);
  
  bool inversePower(Sparse<double>& mat, double eigenEst, std::vector<double>& eigenVector, double& eigenValue, RngWrapper& rng);
  
  bool luDecomp(std::map<std::pair<int,int>, double>& matrix, int nRows, int nCols, std::vector<int>& rowOrder);
  
  bool luSolve(std::map<std::pair<int,int>, double> A, int nRows, int nCols, std::vector<int> rowOrder, std::vector<double>& b, std::vector<double>& X);
  
  bool inversePowerTri(std::vector<double> diag, std::vector<double>& offDiag, double eigenEst, std::vector<double>& eigenVector, double& eigenValue, RngWrapper& rng);
  
  void luSolveTri(const std::vector<double>& diag, const std::vector<double>& lowerOffDiag, const std::vector<double>& upperOffDiag, std::vector<double>& x);
  
  void luDecompTri(std::vector<double>& diag, std::vector<double>& lowerOffDiag, const std::vector<double>& upperOffDiag);
  
  double bisection(const std::vector<double> &alpha, const std::vector<double> &beta, double lower, double upper);
  
  double normL2(const std::vector<double> vect);
  
  void normaliseL2(std::vector<double>& vect);
  
  bool partitionGraph(const Sparse<double>& laplacian, std::vector<int>& groups, int nSplits, int useMedian, RngWrapper& rng);
  std::vector<double> scalerDiv(const std::vector<double>& ,const double);
  std::vector<double> scalerMult(const std::vector<double>& vec, double scaler);
  std::vector<double> mult(const Sparse<double>& mat, const std::vector<double>& vect);
  std::vector<double> mult(const Matrix<double>& mat,const std::vector<double>& vect);
  void scalerSubtract(std::vector<double>& vec, double scaler);
  void scaleVector(std::vector<double>& vec, double scaler);
  
  double dot(const std::vector<double>& a, const std::vector<double>& b);
  std::vector<double> sum(const std::vector<double>& a, const std::vector<double>& b);
  std::vector<double> diff(const std::vector<double>& a, const std::vector<double>& b);
  
  
};


#endif
