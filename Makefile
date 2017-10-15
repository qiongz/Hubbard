CC =  icpc
CFLAGS = -std=c++11 -qopenmp
LIBS= -llapack -lpthread

hubbard:main.cpp basis.o matrix.o init.o lanczos_hamiltonian.o hamiltonian.o
	$(CC) $(CFLAGS) $^ -O2 -o $@ ${LIBS} $(CFLAGS)

basis.o:basis.cpp basis.h
	$(CC) $(CFLAGS) -c basis.cpp

matrix.o:matrix.cpp matrix.h
	$(CC) $(CFLAGS) -c matrix.cpp -o $@

init.o:init.cpp init.h
	$(CC) $(CFLAGS) -c init.cpp

hamiltonian.o:hamiltonian.cpp hamiltonian.h matrix.h
	$(CC) $(CFLAGS) -c hamiltonian.cpp

lanczos_hamiltonian.o:lanczos_hamiltonian.cpp lanczos_hamiltonian.h matrix.h
	$(CC) $(CFLAGS) -c lanczos_hamiltonian.cpp

.PHONY: all clean remove
all: clean hubbard

clean:
	rm -f  *.o 
remove:
	rm hubbard
