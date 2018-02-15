SRC = ./src/

.PHONY: all clean

all:
	@echo "No target is set."

pp:
	make -C $(SRC) $@

clean:
	rm -f *~
	make -C $(SRC) $@
