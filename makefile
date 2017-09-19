CC= g++
CFLAGS=  -Wall -std=c++11 -O3 -g -c -I/usr/local/include
NAME1= sudakov_g
OBJCS = main.o sudakov_g.o Blad.o
NAME= Sudakov_g.exe
all: main

main: $(OBJCS)
	$(CC) -L/usr/local/lib $(OBJCS) -lgsl -lgslcblas -lm -o main.out
main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp
sudakov_g.o: sudakov_g.cpp
	$(CC) $(CFLAGS) sudakov_g.cpp
Blad.o: Blad.cpp
	$(CC) $(CFLAGS) Blad.cpp
clean:
	rm *.o
