all:
	g++ -Wall -shared -fPIC -o solverscpp solverscpp.cpp

clean:
	rm solverscpp
