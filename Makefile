CC=gcc
CFLAGS=-O3
OBJECTS=ASM1.o gpat.o
LINKER_FLAGS=-lm -lpthread

%.o : %.c
	$(CC) $(CFLAGS) -c $^

gpat: $(OBJECTS)
	$(CC) -o $@ $^ $(LINKER_FLAGS)

clean:
	rm -rf gpat *.o
