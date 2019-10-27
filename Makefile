SCRIPT := main.h
EXE := main
CXX := g++
CXXFLAGS := -std=c++11 -O3 -fopenmp
CXX_WARNING_OPTIONS := -Wall -Wextra -Wno-unused-result
LDLIBS := -pthread -lpthread
LIBS := -I ~/lib/eigen/ -I ./spectra/include/ -I ./lib

main: main.cpp matrix.h
	$(CXX) $(LIBS) $(CXXFLAGS) -o main main.cpp matrix.h

clean:
	rm ./$(EXE)
