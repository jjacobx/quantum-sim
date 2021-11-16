CXX = g++
CXXFLAGS = -g -O0 -Wall -Wextra -std=c++20 -mavx2
LDFLAGS = $(shell pkg-config --libs gtest eigen3)
INCS = -I../../src $(shell pkg-config --cflags gtest eigen3)

ODIR = ../../src/obj
ODIR_TEST = obj
_OBJ = circuit.o 
_OBJ_TEST = main_testall.o circuit_test.o
OBJ = $(patsubst %,$(ODIR_TEST)/%,$(_OBJ_TEST)) $(patsubst %,$(ODIR)/%,$(_OBJ)) 

BDIR = ../bin


$(ODIR_TEST)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCS) -c -o $@ $< 

test_all.out: $(OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCS) -o $(BDIR)/$@ $^ 

.PHONY: clean

clean:
	rm -f core $(ODIR_TEST)/*.o

$(shell mkdir -p $(ODIR_TEST))
