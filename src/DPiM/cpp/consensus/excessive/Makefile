CXX = g++ -std=gnu++0x
CXXFLAGS = -Wall 
#CXXFLAGS = -Wall -O3
#CXXFLAGS = -Wall -Wno-deprecated -O3

LIBS = Network.o
GSL = -lgsl -lgslcblas
JSON = -I/home/glocke/DPiM/cpp/rapidjson/rapidjson-master/include
BOOST = -lboost_program_options 

all: ${LIBS} parseNetworks test

Network.o: Network.hpp Network.cpp
	$(CXX) -c Network.cpp $(CXXFLAGS) $(JSON)

parseNetworks: parseNetworks.cpp $(LIBS)
	$(CXX) parseNetworks.cpp -o parseNetworks $(LIBS) $(CXXFLAGS) $(JSON)

test: test.cpp $(LIBS)
	$(CXX) test.cpp -o test $(LIBS) $(CXXFLAGS) $(JSON)

clean:
	rm *.o parseNetworks test a.out
