CXX = g++
CXXFLAGS = -g -Og -Wall -Wextra -std=c++20 -mavx2
LDFLAGS = $(shell pkg-config --libs eigen3)
INCS = $(shell pkg-config --cflags eigen3)

DEPS = circuit.h

ODIR = obj
_OBJ = main.o circuit.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

BDIR = ../bin


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCS) -c -o $@ $< 

sim.out: $(OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCS) -o $(BDIR)/$@ $^

.PHONY: clean

clean:
	rm -f core $(ODIR)/*.o

$(shell mkdir -p $(ODIR))
