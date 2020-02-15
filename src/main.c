#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "utils.h"
#include "genetics.h"
#include "RKF78.h"

/*
 * The following code solves the Gompartz growth model from [McC05, Example 5.3]
 * of finding a minimal solution for a cancer treatment's model.
 * This solution is found using the well-known Gennetic's Algorithm. In this
 * case, the code uses a onPointCrossover method as Crossover function, the
 * BitFlipMutation as mutation method and Tournament Selection with d=2 as
 * Selection With Replacement.
 * The code recivies as inputs the following arguments in the command line:
 *  int popsize : size of the population of individuals
 *  double probmut : probability of mutation in the bitflip method
 *  int typesol (optional) : type of solution method, by default uses 0.
 *      - 0 : the program tries to find a curative solution, if it finds it,
 *      then the code finish, otherwise, it tries to find a paliative solution
 *      - 1 : it tries to find straightly a paliative solution 
 *      - 2 : it will try to find a curative solution, then, whateverit was the
 *      result, it will try to find a paliative solution
 */

int main(int argc, char* argv[]) {

    int popsize, typesol;
    double probmut;
    typesol = 0;

#define NARGS 3
    if (argc<3
        || sscanf(argv[1], "%d", &popsize)!=1
        || sscanf(argv[2], "%lf", &probmut)!=1
       ) {
        fprintf(stderr,"%s popsize probmut [ typesol ]\n", argv[0]);
        return -1;
    }
    if (argc>=NARGS + 1) {
        if (sscanf(argv[NARGS], "%d", &typesol)!=1) {
            fprintf(stderr, "GA(): error while reading type of solution\n");
            return -1;
        }
    }
    if (typesol == 0) {
        /* trying to find a curative solution */
        fprintf(stdout, "# GA() : finding a Curative Solution ... \n");
        if (!geneticAlgorithm(popsize, probmut, Curative_Fitness, writeCurativeSolution)) {
            fprintf(stdout, "# GA() : finding a Paliative Solution ... \n");
            geneticAlgorithm(popsize, probmut, Paliative_Fitness, writePaliativeSolution);
        }
    }
    else if (typesol == 1) {
        fprintf(stdout, "# GA() : finding a Paliative Solution ... \n");
        geneticAlgorithm(popsize, probmut, Paliative_Fitness, writePaliativeSolution);
    }
    else if (typesol == 2) {
        fprintf(stdout, "# GA() : finding a Curative Solution ... \n");
        geneticAlgorithm(popsize, probmut, Curative_Fitness, writeCurativeSolution);
        fprintf(stdout, "# GA() : finding a Paliative Solution ... \n");
        geneticAlgorithm(popsize, probmut, Paliative_Fitness, writePaliativeSolution);
    }
    else {
        fprintf(stderr, "# GA() : type of solution not valid, use 0, 1 or 2\n");
        return -1;
    }
    return 0;
}
