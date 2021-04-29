CXX=c++
CFLAGS = --std=c++17 -O3 -Wall -Wpedantic
# note: if using clang, replace lstdc++fs with lc++fs (libstdc++ vs libc++)
LDFLAGS = -lpng -lstdc++fs -lfmt -fopenmp -lpthread

default: wolff.cpp wolff.hpp main.o
	${CXX} ${CFLAGS} wolff.cpp main.o ${LDFLAGS} -o wolff

main.o: main.cpp wolff.cpp wolff.hpp
	${CXX} ${CFLAGS} -c main.cpp

clean:
	rm *.o wolff
