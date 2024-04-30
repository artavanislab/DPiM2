#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
//#include <cmath> // for sqrt

// adapted from http://www.johndcook.com/blog/standard_deviation/
// "[it] goes back to a 1962 paper by B. P. Welford and is presented in Donald Knuth’s Art of Computer Programming, Vol 2, page 232, 3rd edition."

#ifndef SYM_MAT
#define SYM_MAT 1
typedef boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper> symMat;
#endif

class MatrixRunningStat {
public:
  MatrixRunningStat(size_t size_);
  void Resize(size_t newSize);  
  void Clear();

  void Push(symMat x);

  unsigned NumDataValues() const;
  size_t Size() const;

  symMat Mean() const;
  symMat Variance() const;

private:
  unsigned cnt;
  size_t size;
  symMat prevMean, newMean, prevSqr, newSqr;
};
