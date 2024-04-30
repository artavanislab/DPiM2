CXX = g++ -std=gnu++0x
CXXFLAGS = -Wall -O3
#CXXFLAGS = -Wall -Wno-deprecated -O3

LIBS = Histogram.o BaitByPrey.o PreyByPrey.o 
GSL = -lgsl -lgslcblas
JSON = -I/home/glocke/DPiM/cpp/rapidjson/rapidjson-master/include
BOOST = -L/camhpc/pkg/boost/1.61.0/centos6/lib -lboost_program_options 

all: ${LIBS} hyperspec run_prob

Histogram.o: Histogram.hpp Histogram.cpp
	$(CXX) -c Histogram.cpp $(CXXFLAGS) $(GSL)

BaitByPrey.o: BaitByPrey.hpp BaitByPrey.cpp
	$(CXX) -c BaitByPrey.cpp $(CXXFLAGS) $(GSL)

PreyByPrey.o: PreyByPrey.hpp PreyByPrey.cpp
	$(CXX) -c PreyByPrey.cpp $(CXXFLAGS) $(GSL) $(JSON)

hyperspec: HyperSpec.cpp $(LIBS)
	$(CXX) HyperSpec.cpp -o hyperspec $(LIBS) $(CXXFLAGS) $(GSL) $(BOOST)

run_prob: run_prob.cpp Combinatorics.hpp
	$(CXX) run_prob.cpp -o run_prob $(CXXFLAGS)

clean:
	rm *.o hyperspec a.out
