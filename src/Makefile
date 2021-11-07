CC = g++
CFLAGS = -O0

DEPS = circuit.h

ODIR = obj
_OBJ = main.o circuit.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

BDIR = ../bin


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

sim-test: $(OBJ)
	$(CC) -o $(BDIR)/$@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f core $(ODIR)/*.o
