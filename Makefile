CXX = g++
OBJECTS = main.o 

exec:  $(OBJECTS)
	$(CXX) main.cpp  -o exec 
	
clean :
	rm main.o exec