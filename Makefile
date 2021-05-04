CXX=c++
CFLAGS = --std=c++17 -O3 -Wall -Wpedantic -Wfatal-errors
# note: if using clang, replace lstdc++fs with lc++fs (libstdc++ vs libc++)
LDFLAGS = -lpng -lstdc++fs -lfmt -lpthread -lboost_program_options

default: main.o
	${CXX} ${CFLAGS} sim.cpp main.o ${LDFLAGS} -o ising

main.o: main.cpp
	${CXX} ${CFLAGS} -c main.cpp

clean:
	rm *.o ising
