CXX = g++
OBJECTS = main.o 
FLAGS = -std=c++11 -g -O3 -fopenmp -lpthread

exec:
	$(CXX) $(FLAGS) L-BFGS/lbfgs.c main.cpp -o exec 
	
clean :
	rm $(OBJECTS) exec
