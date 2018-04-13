CC =  icpc
#CFLAGS = -qopenmp -std=c++11
CFLAGS = -openmp
LIBS= -llapack -lpthread

check:check.cpp basis.o matrix.o init.o lanczos_hamiltonian.o hamiltonian.o mt19937-64.o  Greens_function.o
	$(CC) $(CFLAGS) $^ -O2 -o $@ ${LIBS} $(CFLAGS) -lgsl

hubbard:main.cpp basis.o matrix.o init.o lanczos_hamiltonian.o hamiltonian.o mt19937-64.o 
	$(CC) $(CFLAGS) $^ -O2 -o $@ ${LIBS} $(CFLAGS) 

basis.o:basis.cpp basis.h
	$(CC) $(CFLAGS) -c basis.cpp

matrix.o:matrix.cpp matrix.h mt19937-64.h
	$(CC) $(CFLAGS) -c matrix.cpp -o $@

init.o:init.cpp init.h
	$(CC) $(CFLAGS) -c init.cpp

hamiltonian.o:hamiltonian.cpp hamiltonian.h matrix.h
	$(CC) $(CFLAGS) -c hamiltonian.cpp

lanczos_hamiltonian.o:lanczos_hamiltonian.cpp lanczos_hamiltonian.h matrix.h
	$(CC) $(CFLAGS) -c lanczos_hamiltonian.cpp

Greens_function.o:Greens_function.h Greens_function.cpp
	$(CC) $(CFLAGS) -c Greens_function.cpp

mt19937-64.o:mt19937-64.c mt19937-64.h
	$(CC) -c mt19937-64.c 	

.PHONY: all clean remove
all: clean hubbard check

clean:
	rm -f  *.o 
remove:
	rm hubbard check
