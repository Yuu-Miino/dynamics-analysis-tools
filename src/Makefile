.SUFFIXES: .cpp .hpp .o

CC     = g++ -Wall
MEMORY = -fsanitize=address -fno-omit-frame-pointer
OBJS   = Dynamics.o ODE.o Event.o
DIST   = ../bin/

INCLUDE = -I/usr/local/include/ -I../include/

pp:pp.o DynamicsPP.o ODE.o Event.o
	$(CC) -o $(DIST)$@ $^ 

.cpp.o:
	$(CC) -c $(INCLUDE) $<

clean:
	rm -f *~ *.o
