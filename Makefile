CC = g++ 
CFLAGS = -std=c++11 
LIBS= -llapack

hubbard:main.cpp basis.o init.o lanczos.o
	$(CC) $(CFLAGS) $^ -o $@ ${LIBS}

basis.o:basis.cpp basis.h
	$(CC) $(CFLAGS) -c basis.cpp

init.o:init.cpp init.h
	$(CC) $(CFLAGS) -c init.cpp

lanczos.o:lanczos.cpp lanczos.h
	$(CC) $(CFLAGS) -c lanczos.cpp

.PHONY: all clean remove
all: clean hubbard

clean:
	rm -f  *.o 
remove:
	rm hubbard