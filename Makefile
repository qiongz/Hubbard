CC = g++ 
CFLAGS = -std=c++11 
LIBS= -llapack

hubbard:main.cpp basis.o init.o lanczos.o
	$(CC) $(CFLAGS) $^ -o $@ ${LIBS}

basis.o:basis.cpp basis.h
	$(CC) $(CFLAGS) -c $^

init.o:init.cpp init.h
	$(CC) $(CFLAGS) -c $^

lanczos.o:lanczos.cpp
	$(CC) $(CFLAGS) -c $^

.PHONY: all clean remove
all: clean hubbard

clean:
	rm -f  *.o 
remove:
	rm hubbard
