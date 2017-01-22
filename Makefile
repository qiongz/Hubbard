CC = g++ 
CFLAGS = -std=c++11 
LIBS= -llapack

hubbard:main.cpp basis.o init.o lanczos.o hamiltonian.o
	$(CC) $(CFLAGS) $^ -o $@ ${LIBS}

basis.o:basis.cpp basis.h
	$(CC) $(CFLAGS) -c basis.cpp

init.o:init.cpp init.h
	$(CC) $(CFLAGS) -c init.cpp

hamiltonian.o:hamiltonian.cpp hamiltonian.h
	$(CC) $(CFLAGS) -c hamiltonian.cpp

lanczos.o:lanczos.cpp lanczos.h
	$(CC) $(CFLAGS) -c lanczos.cpp

hamiltonian.o:hamiltonian.cpp hamiltonian.h
	$(CC) $(CFLAGS) -c hamiltonian.cpp

.PHONY: all clean remove
all: clean hubbard

clean:
	rm -f  *.o 
remove:
	rm hubbard
