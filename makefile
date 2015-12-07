CC=gcc
CFLAGS=-Wall -ansi
LDFLAGS=-lgmp
EXEC=bin/main
SRC=src/main.c src/sylvester.c src/lagrange.c src/horner.c src/tools.c src/thm_chinois.c src/solution.c
OBJ=obj/main.o obj/sylvester.o obj/lagrange.o obj/horner.o obj/tools.o obj/thm_chinois.o obj/solution.o

all: $(EXEC)

bin/main : $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

obj/main.o : $(SRC)
	$(CC) -c src/main.c -o obj/main.o $(CFLAGS)

obj/%.o : src/%.c
	$(CC) -o $@ -c $< $(CFLAGS)

clean :
	rm -f obj/*.o bin/main
