.SUFFIXES: .cpp .hpp .o

CC      = g++ -Wall
LIBS    = -lm
INCLUDE = -I/usr/local/include/ -I../include/
MEMORY  = -fsanitize=address -fno-omit-frame-pointer
DIST    = ../bin/

.PHONY: all clean

all:
	@echo "No target is set."

pp:pp.o DynamicsPP.o ODE.o Event.o
	$(CC) -o $(DIST)$@ $(OBJDIR)$^ 

.cpp.o:
	$(CC) -o $(INCLUDE) $@ -c $<

clean:
	rm -f *~ *.o
