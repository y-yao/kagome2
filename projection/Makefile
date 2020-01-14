EXE := project2
CXX := g++
CXXFLAGS := -std=c++11 -O3 -fopenmp
CXX_WARNING_OPTIONS := -Wall -Wextra -Wno-unused-result
LDLIBS := -pthread -lpthread
LIBS := -I ../lib

$(EXE): dimers2.cpp
	$(CXX) $(LIBS) $(CXXFLAGS) -o $(EXE) dimers2.cpp

clean:
	rm ./$(EXE)
