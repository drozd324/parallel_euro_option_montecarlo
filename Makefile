CC :=  gcc #clang++
CCFLAGS := -lm -lgsl -Wextra -Wall 
DEBUG = $(" ") #-g -fsanitize=address -Wall -Wextra -lefence #$(" ") #

execs = main parallel_main

all: $(execs)

parallel_main: parallel_main.c
	$(CC) -o $@ $^ $(CCFLAGS) $(DEBUG) -fopenmp

main: main.c
	$(CC) -o $@ $^ $(CCFLAGS) $(DEBUG)

.PHONY: clean

clean:
	rm -f *.o $(execs)
