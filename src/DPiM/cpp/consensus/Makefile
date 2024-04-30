CXX = g++ -std=gnu++0x
#CXXFLAGS = -Wall -g -O3
CXXFLAGS = -Wall -O3
#CXXFLAGS = -Wall -Wno-deprecated -O3

LIBS = MatrixRunningStat.o Network.o
GSL = -lgsl -lgslcblas
JSON = -I/home/glocke/DPiM/cpp/rapidjson/rapidjson-master/include
BOOST = -lboost_program_options 

all: ${LIBS} parseNetworks0 dum test

MatrixRunningStat.o: MatrixRunningStat.hpp MatrixRunningStat.cpp
	$(CXX) -c MatrixRunningStat.cpp $(CXXFLAGS) 

Network.o: Network.hpp Network.cpp MatrixRunningStat.o 
	$(CXX) -c Network.cpp  $(CXXFLAGS) $(JSON)

parseNetworks0: parseNetworks0.cpp $(LIBS)
	$(CXX) parseNetworks0.cpp -o parseNetworks0 $(LIBS) $(CXXFLAGS) $(JSON) $(BOOST)

test: test.cpp $(LIBS)
	$(CXX) test.cpp -o test $(LIBS) $(CXXFLAGS) $(JSON)

dum: dum.cpp $(LIBS)
	$(CXX) dum.cpp -o dum $(LIBS) $(CXXFLAGS) $(JSON)

clean:
	rm *.o parseNetworks0 test dum a.out
