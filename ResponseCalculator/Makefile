LIBS = $(shell root-config --libs)
INCS = -I$(shell root-config --incdir) -I.

WARN = -Wall -Wshadow

OBJ = ResponseCalculator.o
EXE = ResponseCalculator

all: $(EXE) Makefile
	@

%.o:%.cc
	g++ -o $@ $+ -c -O2 $(INCS) $(WARN) -std=c++11

$(EXE): $(OBJ)
	g++ -O2 -o $@ $+ $(LIBS) $(WARN)


clean:
	@rm -f $(OBJ)
	@rm -f *.cxx lib*h # remove rootcints generated source and header files
	@rm -f AutoDict_* *_h.d *_cc.d *_C.d # cint generated files


