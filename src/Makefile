.SUFFIXES: .cpp .hpp .o

CC      = g++ -Wall
LIBS    = -lm
INCLUDE = -I/usr/local/include/ -I../include/
MEMORY  = -fsanitize=address -fno-omit-frame-pointer
DIST    = ../bin/
BASE_O  = HSODEsolver.o ODEsolver.o pwlDuffing.o pwlDuffingMode.o 

.PHONY: all clean

all:
	@echo "No target is set."

pp: pp.o $(BASE_O)
	$(CC) -o $(DIST)$@ $(OBJDIR)$^ $(LIBS)
ppPendulum: ppPendulum.o pendulum.o ODEsolver.o
	$(CC) -o $(DIST)$@ $(OBJDIR)$^ $(LIBS)
bf1: bf1.o $(BASE_O)
	$(CC) -o $(DIST)$@ $(OBJDIR)$^ $(LIBS)
fix: fix.o $(BASE_O)
	$(CC) -o $(DIST)$@ $(OBJDIR)$^ $(LIBS)
bf2: bf2.o $(BASE_O)
	$(CC) -o $(DIST)$@ $(OBJDIR)$^ $(LIBS)

.cpp.o:
	$(CC) $(INCLUDE) -c $<

clean:
	rm -f ./*~ ../include/*~ *.o
