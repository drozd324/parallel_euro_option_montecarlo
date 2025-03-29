CC :=  gcc #clang++
CCFLAGS := -lm -Wextra -Wall -Wsign-conversion -Werror --std=c++23
DEBUG = -g -fsanitize=address -Wall -Wextra -lefence #$(" ") #

all: main 

main: main.cc
	$(CC) $(CCFLAGS) -o $@ $^ $(DEBUG)

%.o: %.cc
	$(CC) $(CCFLAGS) -c $< -o $@

.PHONY: clean

clean:
	rm -f *.o $(execs)
