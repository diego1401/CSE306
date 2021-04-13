CXX = g++
OBJECTS = main.o 

exec:  $(OBJECTS)
	$(CXX) -O3 main.cpp  -o exec 
	
clean :
	rm main.o exec