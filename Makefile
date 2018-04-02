SRC = ./src/
BIN = ./bin/
INC = ./include/

.PHONY: all clean
.IGNORE: clean

all:
	@echo "No target is set."

pp:
	make -C $(SRC) $@
bf1:
	make -C $(SRC) $@
fix:
	make -C $(SRC) $@
bf2:
	make -C $(SRC) $@

clean:
	rm -f ./*~
	make -C $(SRC) $@
	make -C $(BIN) $@
	make -C $(INC) $@
