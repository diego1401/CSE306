CXX = g++-10
OBJECTS = main.o 
FLAGS =  -g -O3 -fopenmp -lpthread

exec:
	$(CXX) $(FLAGS) L-BFGS/lbfgs.c main.cpp -o exec 
	
clean :
	rm $(OBJECTS) exec
