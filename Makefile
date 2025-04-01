CC =  gcc #clang++
MPICC = mpicc
CCFLAGS = -O0 -lm -lgsl -Wextra -Wall 
DEBUG = $(" ")# -g -fsanitize=address -Wall -Wextra -lefence #$(" ") #

EXECS = main parallel_main_mpi parallel_main_omp 

all: $(EXECS)

parallel_main_mpi: parallel_main_mpi.c
	$(MPICC) -o $@ $^ $(CCFLAGS) $(DEBUG) 

parallel_main_omp: parallel_main_omp.c
	$(CC) -o $@ $^ $(CCFLAGS) $(DEBUG) -fopenmp

main: main.c
	$(CC) -o $@ $^ $(CCFLAGS) $(DEBUG) 

.PHONY: clean

clean:
	rm -f *.o $(EXECS)
