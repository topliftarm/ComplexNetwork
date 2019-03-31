CXX=g++
CXXFLAGS  = -std=c++11 -O2 -fPIC -c -fopenmp
SWIGFLAGS = -c++ -python -shadow  

# location of the Python header files
PYTHON_VERSION = 3.7
PYTHON_INCLUDE = /home/vahid/Programming/anaconda3/include/python3.7m


TARGET = ode_solver
OBJS= $(TARGET).o $(TARGET)_wrap.o
HFILES= ode_solver.h 

_$(TARGET).so: $(OBJS)
	$(CXX) -shared -fopenmp $(OBJS) -o _$(TARGET).so
 
$(TARGET)_wrap.o: $(TARGET)_wrap.cpp $(HFILES)
	$(CXX) $(CXXFLAGS) $(TARGET)_wrap.cpp -I $(PYTHON_INCLUDE) 

$(TARGET).o: $(TARGET).cpp  $(TARGET).h
	$(CXX) $(CXXFLAGS) $(TARGET).cpp

$(TARGET)_wrap.cpp : $(TARGET).i 
	swig $(SWIGFLAGS) -o $(TARGET)_wrap.cpp $(TARGET).i 

clean :
	rm  ../data/f/*

eradicate:
	rm *.o  *.so $(TARGET)_wrap.cpp $(TARGET).py *.pyc ../data/f/* 

.PHONEY : clean eradicate




# swig -c++ -python  -shadow number.i
# g++ -O2 -fPIC -c number.cxx
# g++ -O2 -fPIC -c number_wrap.cxx -I /usr/include/python2.7 
# g++ -shared number.o number_wrap.o -o _number.so
