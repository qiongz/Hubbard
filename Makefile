CC = g++ 
CFLAGS = -std=c++11

hubbard:main.cpp basis.o
	$(CC) $(CFLAGS) $^ -o $@

basis.o:basis.cpp basis.h
	$(CC) $(CFLAGS) -c $^

.PHONY: all clean remove
all: clean hubbard

clean:
	rm -f  *.o 
remove:
	rm hubbard
