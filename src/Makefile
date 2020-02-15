# Optimization
OPT=-g -Wall

# All
all: main
superall: main solutions

solutions: main
	./main 40=popsize 0.1=probmutation 2=typesol
	@echo 1>&2 "*** Application of the Genetic Algorithm"
	@echo 1>&2 "*** Try in gnuplot:"
	@echo 1>&2 "*** plot 'curative_sol.sol' u 1:2 w l"

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
	rm -f main *.sol
