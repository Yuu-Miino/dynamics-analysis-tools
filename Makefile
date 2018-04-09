SRC = ./src/
BIN = ./bin/
INC = ./include/

.PHONY: all clean
.IGNORE: clean

all: pp ppPendulum bf1 fix bf2 gz2

pp:
	make -C $(SRC) $@
ppPendulum:
	make -C $(SRC) $@
bf1:
	make -C $(SRC) $@
fix:
	make -C $(SRC) $@
bf2:
	make -C $(SRC) $@
gz2:
	make -C $(SRC) $@

clean:
	rm -f ./*~
	make -C $(SRC) $@
	make -C $(BIN) $@
	make -C $(INC) $@
