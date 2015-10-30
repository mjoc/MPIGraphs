#ifndef __MOC_SparseBinary_
#define __MOC_SparseBinary_
#include <vector>
class Runlog;


class SparseBinary{
public:
  int _nRows;
  int _nCols;
  std::vector<int> _cols;
  std::vector<int> _rowIndex;
  SparseBinary();
	SparseBinary(int nRows, int nCols, const std::vector<int>& colsArray, const std::vector<int>& rowStartIdxArray);
  SparseBinary(const SparseBinary&);
  SparseBinary& operator=(const SparseBinary&);
  void print();
  void print(Runlog& runlog);
  bool readFromFile(char *filename);
  int nRows()const {return _nRows;};
  int nCols()const {return _nCols;};
  int nData()const {return _rowIndex[_nRows];}
  int operator()(int iRow, int iCol);
};




#endif
