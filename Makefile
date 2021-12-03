CXX=g++
CXXFLAGS=-Wall -Wextra -pedantic -std=c++17 -O0 -g
LDFLAGS=$(CXXFLAGS)
OBJ=$(SRC:.cc=.o)

all:  tsp

tsp: tsp.o chromosome.o deme.o cities.o
	$(CXX) $(LDFLAGS) -o $@ $^

deme_test: deme_test.o chromosome.o deme.o cities.o
	$(CXX) $(LDFLAGS) -o $@ $^

%.o: %.cc %.hh
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) -c -o $@ $<

clean:
	rm -rf *.o *.out *.gch tsp deme_test chromosome_test
