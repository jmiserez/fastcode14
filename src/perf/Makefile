CC     = gcc
CFLAGS = -O3 -g -std=c99
LIBS   = -lpfm
SRC    = example.c other_comp.c

example_with_counters: $(SRC)
	$(CC) $(CFLAGS) $(LIBS) -DCOST_MEASURE -o $@ $^

example_without_counters: $(SRC)
	$(CC) $(CFLAGS) $(LIBS) -o $@ $^

all: example_with_counters example_without_counters

run:
	@echo
	@echo === Example with counters disabled ==============================
	@echo
	@./example_without_counters
	@echo
	@echo === Example with counters enabled ===============================
	@echo
	@./example_with_counters
