# Dir setting
HOME			= .
INCLUDE 		= $(HOME)/include
SRC 			= $(HOME)/src
BIN				= $(HOME)/bin
LIB				= $(HOME)/lib


# Target & target specific main file
TARGET			= $(BIN)/Shot
MAIN_TAR		= $(SRC)/main.cpp

# Sources
SOURCES			= $(wildcard $(SRC)/*.cpp) $(wildcard $(SRC)/*/*.cpp)
OBJECTS			= $(patsubst %.cpp, %.o, $(SOURCES))

# Compiler
RM				= rm -rf
CXX				= g++
CXXFLAGS		= -std=c++11 -g -fmessage-length=0 -c -fopenmp -msse3 -mavx -I${INCLUDE}/ARPACKpp -Ofast -ffast-math
LDFLAGS			= -fopenmp
DEP				= $(LIB)/libarpack.a $(LIB)/libopenblas.a
LDLIBS			= -lgfortran -lm 


## Binaries
$(TARGET): $(OBJECTS) $(DEP)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.cpp $(LIB)/libarpack.a
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(LIB)/libopenblas.a:
	cd $(LIB) && sh install-openblas.sh

$(LIB)/libarpack.a: $(LIB)/libopenblas.a
	cd $(LIB) && sh install-arpack-ng.sh

# Public rules
all: $(TARGET)

clean:
	$(RM) $(OBJECTS) $(TARGET)

clear:
	$(RM) $(OBJECTS)

.PHONY: all clean clear
