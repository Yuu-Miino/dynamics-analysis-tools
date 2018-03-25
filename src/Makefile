.SUFFIXES: .cpp .hpp .o

CC      = g++ -Wall
LIBS    = -lm
INCLUDE = -I/usr/local/include/ -I../include/
MEMORY  = -fsanitize=address -fno-omit-frame-pointer
DIST    = ../bin/

.PHONY: all clean

all:
	@echo "No target is set."

pp: pp.o HSODEsolver.o ODEsolver.o pwlDuffing.o pwlDuffingMode.o
	$(CC) -o $(DIST)$@ $(OBJDIR)$^ $(LIBS)

.cpp.o:
	$(CC) $(INCLUDE) -c $<

clean:
	rm -f ./*~ ../include/*~ *.o
