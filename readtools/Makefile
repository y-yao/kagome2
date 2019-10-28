EXE := read
CXX := g++
CXXFLAGS := -std=c++11 -O3 -fopenmp
CXX_WARNING_OPTIONS := -Wall -Wextra -Wno-unused-result
LDLIBS := -pthread -lpthread
LIBS := -I ../lib

$(EXE): read.cpp
	$(CXX) $(LIBS) $(CXXFLAGS) -o $(EXE) read.cpp

clean:
	rm ./$(EXE)
