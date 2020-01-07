EXE := make_hh
CXX := g++
CXXFLAGS := -std=c++11 -O2 -fopenmp
CXX_WARNING_OPTIONS := -Wall -Wextra -Wno-unused-result
LDLIBS := -pthread -lpthread
LIBS := -I ../lib

$(EXE): make_hh.cpp
	$(CXX) $(LIBS) $(CXXFLAGS) -o $(EXE) make_hh.cpp

clean:
	rm ./$(EXE)
