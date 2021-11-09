CXX = g++
CXXFLAGS = -O0 -Wall -Wextra -std=c++20

DEPS = circuit.h

ODIR = obj
_OBJ = main.o circuit.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

BDIR = ../bin


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

sim.out: $(OBJ)
	$(CXX) -o $(BDIR)/$@ $^ $(CXXFLAGS)

.PHONY: clean

clean:
	rm -f core $(ODIR)/*.o
