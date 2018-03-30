SRC = ./src/

.PHONY: all clean

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
