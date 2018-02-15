SRC = ./src/

all:

pp:
	make -C $(SRC) $@

clean:
	rm -f *~
	make -C $(SRC) $@
