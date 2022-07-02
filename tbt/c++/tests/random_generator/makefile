CXXFLAGS = -std=c++11 -g -Wall --pedantic

distributions: distributions.cpp
	g++ distributions.cpp $(CXXFLAGS) -o distributions.x

plot distributions: 
	./distributions.x 100
	000

