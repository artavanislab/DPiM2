#include "MatrixRunningStat.hpp"

namespace blas = boost::numeric::ublas;

MatrixRunningStat::MatrixRunningStat(size_t setSize) {
  cnt = 0;
  size = setSize;
  prevMean = symMat(size, size);
  newMean = symMat(size, size);
  prevSqr = symMat(size, size);
  newSqr = symMat(size, size);
}

void MatrixRunningStat::Resize(size_t newSize) {
  size = newSize;
  prevMean.resize(size);
  newMean.resize(size);
  prevSqr.resize(size);
  newSqr.resize(size);
}

// doesn't do the work of clearing things, just makes it so that when new data are
//   added, the old data is forgotten.  (that is those memories are cleared only then.)
void MatrixRunningStat::Clear() {
  cnt = 0;
}

// add new data, updating statistics
void MatrixRunningStat::Push(symMat x) {
  cnt++;

  // See Knuth TAOCP vol 2, 3rd edition, page 232
  if (cnt == 1) {
    prevMean = newMean = x;
    prevSqr.clear();
  } else {
    newMean = prevMean + (x - prevMean)/cnt;
    newSqr = prevSqr + blas::element_prod(x - prevMean, x - newMean);
    
    // set up for next iteration
    prevMean = newMean; 
    prevSqr = newSqr;
  }
}

unsigned MatrixRunningStat::NumDataValues() const {
  return cnt;
}
size_t MatrixRunningStat::Size() const {
  return size;
}


symMat MatrixRunningStat::Mean() const {
  return (cnt > 0) ? newMean : 0.0;
}

symMat MatrixRunningStat::Variance() const {
  assert(cnt > 1);
  return newSqr/(cnt - 1);
}
