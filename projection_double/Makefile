EXE := project_double
CXX := g++
CXXFLAGS := -std=c++11 -O3 -fopenmp
CXX_WARNING_OPTIONS := -Wall -Wextra -Wno-unused-result
LDLIBS := -pthread -lpthread
LIBS := -I ../lib

$(EXE): dimers_double.cpp
	$(CXX) $(LIBS) $(CXXFLAGS) -o $(EXE) dimers_double.cpp

clean:
	rm ./$(EXE)
