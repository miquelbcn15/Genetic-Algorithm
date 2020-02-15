# Optimization
OPT=-g -Wall

# All
all: main

main: main.c genetics.o utils.o RKF78.o
	gcc -o main $(OPT) main.c genetics.o utils.o RKF78.o -lm

# Libraries
genetics.o: genetics.c genetics.h 
	gcc -c  $(OPT) genetics.c 
utils.o: utils.c utils.h
	gcc -c $(OPT) utils.c
RKF78.o: RKF78.c RKF78.h
	gcc -c $(OPT) RKF78.c

# Clean
clean: 
	rm -f *.o
realclean: clean
	rm -f main 
