#include "MatrixRunningStat.hpp"

namespace blas = boost::numeric::ublas;

MatrixRunningStat::MatrixRunningStat() {
  init(0);
}
MatrixRunningStat::MatrixRunningStat(size_t setSize) {
  init(setSize);
}
void MatrixRunningStat::init(size_t setSize) {
  cnt = 0;
  size_ = setSize;
  varFresh=false;

  prevMean = symMat(size_, size_, size_*(size_-1)/2);
  newMean = symMat(size_, size_, size_*(size_-1)/2);
  newSqr = symMat(size_, size_, size_*(size_-1)/2);
  prevSqr = symMat(size_, size_, size_*(size_-1)/2);
}

void MatrixRunningStat::resize(size_t newSize) {
  std::cerr << "MatrixRunningStat::resize is not implemented fully. "
	    << "It must be changed to zero-out the new slices." << std::endl;
  size_ = newSize;
  prevMean.resize(size_, size_); 
  newMean.resize(size_, size_);
  prevSqr.resize(size_, size_);
  newSqr.resize(size_, size_);
}

void MatrixRunningStat::clear() {
  cnt = 0;
  prevMean.clear();
  newMean.clear();
  prevSqr.clear();
  newSqr.clear();
}

// add new data, updating statistics
void MatrixRunningStat::push(symMat x) {
  assert(x.size1() == size_);
  cnt++;
  varFresh=false;
  
  
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

unsigned MatrixRunningStat::getCnt() const {
  return cnt;
}

size_t MatrixRunningStat::size() const {
  return size_;
}

symMat MatrixRunningStat::mean() const {
  assert(cnt > 0);
  return newMean;
}
double MatrixRunningStat::mean(size_t row, size_t col) const {
  assert(cnt > 0);
  return newMean(row, col);
}

symMat MatrixRunningStat::variance() {
  if (!varFresh) {
    assert(cnt > 1);
    calcVariance();
  }
  return var;
}
double MatrixRunningStat::variance(size_t row, size_t col) {
  if (!varFresh) {
    assert(cnt > 1);
    calcVariance();
  }
  return var(row, col);
}

void MatrixRunningStat::calcVariance() {
  if (varFresh) {
    return;
  }
  var = newSqr/(cnt - 1);
  varFresh = true;
}

void MatrixRunningStat::dump() {
  std::cout << "cnt = " << cnt << std::endl;
  std::cout << "newMean = " << newMean << std::endl;
  std::cout << "newSqr = " << newSqr << std::endl;
}
