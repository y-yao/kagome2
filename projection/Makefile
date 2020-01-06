EXE := project
CXX := g++
CXXFLAGS := -std=c++11 -g -fopenmp
CXX_WARNING_OPTIONS := -Wall -Wextra -Wno-unused-result
LDLIBS := -pthread -lpthread
LIBS := -I ../lib

$(EXE): dimers.cpp
	$(CXX) $(LIBS) $(CXXFLAGS) -o $(EXE) dimers.cpp

clean:
	rm ./$(EXE)
