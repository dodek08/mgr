CC= g++
CFLAGS=  -Wall -std=c++11 -O3 -g -c -I/net/software/local/lhapdf/6.2.1/include -I/net/software/local/gsl/2.4/include
NAME1= sudakov_g
OBJCS = main.o sudakov_g.o Blad.o sudakov_updf.o sudakov_f2.o sudakov_cs.o
NAME= Sudakov_g.exe
all: main

main: $(OBJCS)
	$(CC) $(OBJCS)  -L/net/software/local/lhapdf/6.2.1/lib -lLHAPDF -L/net/software/local/gsl/2.4/lib -lgsl -lgslcblas -lm -o main.out
main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp
sudakov_updf.o: sudakov_updf.cpp
	$(CC) $(CFLAGS) sudakov_updf.cpp
sudakov_g.o: sudakov_g.cpp
	$(CC) $(CFLAGS) sudakov_g.cpp
sudakov_f2.o: sudakov_f2.cpp
	$(CC) $(CFLAGS) sudakov_f2.cpp
sudakov_cs.o: sudakov_cs.cpp
	$(CC) $(CFLAGS) sudakov_cs.cpp
Blad.o: Blad.cpp
	$(CC) $(CFLAGS) Blad.cpp
clean:
	rm *.o