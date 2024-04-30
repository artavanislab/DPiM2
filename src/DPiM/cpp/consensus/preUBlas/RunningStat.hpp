#include <cmath> // for sqrt
// taken from http://www.johndcook.com/blog/standard_deviation/


class RunningStat {
public:
  RunningStat() : m_n(0) {}

  // initialize with n elements, only one of which is non-zero
  RunningStat(int n, double x) {
    // this implementation is assured of working, but in principle
    // it would be dead simple to work out how to simply initialize everything
    m_n = 0;
    Push(x);
    while (m_n < n) { Push(0); }
  }
  
  void Clear() {
    m_n = 0;
  }

  void Push(double x) {
    m_n++;

    // See Knuth TAOCP vol 2, 3rd edition, page 232
    if (m_n == 1) {
	m_oldM = m_newM = x;
	m_oldS = 0.0;
    } else {
	m_newM = m_oldM + (x - m_oldM)/m_n;
	m_newS = m_oldS + (x - m_oldM)*(x - m_newM);
    
	// set up for next iteration
	m_oldM = m_newM; 
	m_oldS = m_newS;
      }
  }

  int NumDataValues() const {
    return m_n;
  }

  double Mean() const {
    return (m_n > 0) ? m_newM : 0.0;
  }

  double Variance() const {
    return ( (m_n > 1) ? m_newS/(m_n - 1) : 0.0 );
  }

  double StandardDeviation() const {
    return sqrt( Variance() );
  }

private:
  int m_n;
  double m_oldM, m_newM, m_oldS, m_newS;
};
