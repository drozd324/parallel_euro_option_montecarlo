CC =  gcc #clang++
CCFLAGS = -O3 -lm -lgsl -Wextra -Wall 
DEBUG = $(" ")# -g -fsanitize=address -Wall -Wextra -lefence #$(" ") #

EXECS = main parallel_main

all: $(EXECS)

parallel_main: parallel_main.c
	$(CC) -o $@ $^ $(CCFLAGS) $(DEBUG) -fopenmp

main: main.c
	$(CC) -o $@ $^ $(CCFLAGS) $(DEBUG)

.PHONY: clean

clean:
	rm -f *.o $(EXECS)
