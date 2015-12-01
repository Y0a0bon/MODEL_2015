CC=gcc
CFLAGS=-Wall -ansi
LDFLAGS=-lgmp
EXEC=bin/main
SRC=$(src/)
OBJ=$(obj/)

all: $(EXEC)

bin/main : obj/main.o obj/sylvester.o obj/lagrange.o obj/horner.o
	$(CC) obj/main.o obj/sylvester.o obj/lagrange.o obj/horner.o -o bin/main $(LDFLAGS)

obj/main.o : src/main.c src/sylvester.c src/horner.c
	$(CC) -c src/main.c -o obj/main.o $(CFLAGS)

obj/sylvester.o : src/sylvester.c
	$(CC) -c src/sylvester.c -o obj/sylvester.o $(CFLAGS)

obj/lagrange.o : src/lagrange.c
	$(CC) -c src/lagrange.c -o obj/lagrange.o $(CFLAGS)

obj/horner.o : src/horner.c
	$(CC) -c src/horner.c -o obj/horner.o $(FLAGS)

clean :
	rm -f obj/*.o bin/main
