CXX = g++-10
OBJECTS = main.o 
FLAGS = -g -O3 -fopenmp -lpthread 

exec:  $(OBJECTS)
	$(CXX) $(FLAGS) main.cpp -o exec 
	
clean :
	rm $(OBJECTS) exec
