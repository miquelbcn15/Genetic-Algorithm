#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <float.h>
#include "utils.h"
#include "RKF78.h"
#include "genetics.h"

double uniform(void){ 
    /* function from stackoverflow, uses urandom file which stores entropy of the machine */
    int randomData;
    if ( (randomData = open("/dev/urandom", O_RDONLY)) < 0 ) {
        ExitError("Error while openning urandom file", 1);
    }
    unsigned char randomNumber;
    ssize_t result;
    if ( (result = read(randomData, &randomNumber, 1)) < 0 ){
        close(randomData);
        ExitError("Error while extracting random number", 1);
    }
    close(randomData);
    return randomNumber/255.0;
}

double ran1(long *idum) {
    /* random number generator f Park and Miller from Numerical Recipies */
    /* Call with idum a negative integer to initialize */
    int j;
    long k;
    static long iy = 0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0 || !iy) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        for (j = NTAB+7; j >= 0; j--) {
            k=(*idum)/IQ;
            *idum=IA*(*idum-k*IQ)-IR*k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    if (*idum < 0) *idum += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = *idum;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

void initializePopulation(individual** population, int popsize, fitness_function f) {
    /* Allocating memory for initial population */
    register unsigned char i,k;
    if ( ((*population) = (individual*)malloc(popsize * sizeof(individual))) == NULL) {
        ExitError("when allocating initial population", 23);
    }
    /* Generating random numbers and fill the population */
    unsigned char *Cij;
    if ( (Cij = (unsigned char*)malloc(d_par*(n_par - 1)*popsize*sizeof(unsigned char))) == NULL) {
        ExitError("when allocating super random vector", 23);
    }
    int randomData;
    ssize_t result;
    if ( (randomData = open("/dev/urandom", O_RDONLY)) < 0) {
        ExitError("when opening urandom file", 1);
    }
    if ( (result = read(randomData, Cij, (n_par - 1)*d_par * popsize)) != (n_par - 1)*d_par*popsize ) {
        ExitError("when generating super random vector", 1);
    }
    for ( k = 0; k < popsize; k++) {
        for( i = 0; i < (n_par  - 1)*d_par; i++) (*population)[k].Cij[i] = *(Cij++) >> 4;
        /* Compute fitness function for each individual of the population */
        (*population)[k].fitness = f((*population)[i].Cij);
    }
    close(randomData);
}

void mutation(individual *ind, double probmutation) {
    /* Mutation method based on BitFlipMutation and extracted from */
    /* 18.2.1 Bitwise Operators of K&R Sec. 2.9 */
    for (int i = 0; i < (n_par - 1) * d_par; i++) {
        if (uniform() < probmutation) {
            unsigned char position = uniform() * (4 - 1);
            ind->Cij[i] = ind->Cij[i] ^ (1U << position);
        }
    }
}

void OnePointCrossover(individual* indivp1, individual* indivp2,
                       individual* indivf1, individual* indivf2) {
    /* Crossover method selected is OnePointCrossover from */
    /* 18.2.1 BitWise Operators of K&R Sec. 2.9 */
    for (int i = 0; i < (n_par - 1)* d_par; i++) {
        unsigned char d = uniform() * 4 - 1;
        unsigned char mask = 0xFFFFFFFFU << d;
        indivf1->Cij[i] = (indivp1->Cij[i] & mask) | (indivp2->Cij[i] & ~mask);
        indivf2->Cij[i] = (indivp2->Cij[i] & mask) | (indivp1->Cij[i] & ~mask);
    }
}

individual* TournamentSelection(individual **population, int popsize) {
    /* Tournament selection with t=2, based on the book of Essentials*/
    /* Metaheuristics */
    int firstindex;
    int secondindex;
    do {
        firstindex  = uniform() * (popsize - 1);
        secondindex = uniform() * (popsize - 1);
    } while (firstindex == secondindex);
    if ((*population)[firstindex].fitness < (*population)[secondindex].fitness) {
        return (*population + firstindex);
    }
    return (*population + secondindex);
}


void findFittest(individual* fittest, individual* population, int popsize) {
    double minimum = MAXDOUBLE;
    int index = 0;
    for(int k = 0; k < popsize; k++) {
        if ( minimum > population[k].fitness) index = k;
    }
    fittest->fitness = population[index].fitness;
    for (int i = 0; i < (n_par - 1)*d_par; i++) fittest->Cij[i] = population[index].Cij[i];
}

int geneticAlgorithm(int popsize, double probmutation, fitness_function f, write_function writesol) {
    if ( popsize < 2 
            || !probmutation
            || (popsize % 2)!=0
            ) {
        fprintf(stderr, "Genetic Algorithm cannot work\n");
        return -1;
    }

    /* Initialize initial population of size: popsize and its fitness */
    individual *population;
    initializePopulation(&population, popsize, f);

    /* Finding the best individual inside population: the fittest one Globally */
    individual fittestGlobal;
    findFittest(&fittestGlobal, population, popsize);
    /* Updating the current fittest individual */
    individual fittestActual;
    fittestActual.fitness = fittestGlobal.fitness;
    for (int i = 0; i < (n_par - 1)*d_par; i++) fittestActual.Cij[i] = fittestGlobal.Cij[i]; 

    /* Initializing new Generation */
    individual* nextGeneration;
    if ( (nextGeneration = (individual*)malloc(popsize * sizeof(individual))) == NULL) {
        ExitError("when allocating Next Generation", 23);
    }

    int itermax, iter, globaliter;
    if (fittestActual.fitness == MAXDOUBLE) itermax = 10;
    else itermax = 4;
    iter = 0; globaliter = 0;
    clock_t start_t = clock();
    
    while (iter < itermax) {
        /* Main loop of the Algorithm */
        for (int i = 0; i < popsize/2; i++) {
            /* selection with replacement */
            individual* firstparent  = TournamentSelection(&population, popsize);
            individual* secondparent = TournamentSelection(&population, popsize);
            /* performing crossover */
            OnePointCrossover(firstparent, secondparent, nextGeneration + 2*i,
                              nextGeneration + 2*i + 1);
            /* mutating */
            mutation(nextGeneration + 2*i, probmutation); 
            mutation(nextGeneration + 2*i + 1, probmutation);
            /* Computing their Fitness function */
            (nextGeneration + 2*i)->fitness = f((nextGeneration + 2*i)->Cij);
            (nextGeneration + 2*i + 1)->fitness = f((nextGeneration + 2*i + 1)->Cij);
        }
        /* fittest individual in actual population */
        findFittest(&fittestActual, nextGeneration, popsize);
        if (fittestActual.fitness < fittestGlobal.fitness) {
            fittestGlobal.fitness = fittestActual.fitness;
            for (int i = 0; i < (n_par - 1)*d_par; i++) fittestGlobal.Cij[i] = fittestActual.Cij[i];
            itermax = 5; iter = 0;
        }
        /* Updating new population and cleaning memory */
        individual* aux;
        aux = population;
        population = nextGeneration;
        nextGeneration = aux;
        iter++; globaliter++;
    }
    clock_t end_t = clock();
    free(population); free(nextGeneration);
    if (fittestGlobal.fitness < MAXDOUBLE) {
        fprintf(stdout, "# GA() : solution found in %d iterations and %lf seconds\n",
                globaliter, (double)(end_t - start_t)/CLOCKS_PER_SEC);
        fprintf(stdout, "# GA() : fittest value %lf\n", fittestGlobal.fitness);
        fprintf(stdout, "# GA() : printing solution ... \n");
        writesol(fittestGlobal.Cij);
        return 1;
    }
    else {
        fprintf(stdout, "# GA() : solution not found in %d iterations and %lf seconds\n", 
                globaliter, (double)(end_t - start_t)/CLOCKS_PER_SEC);
        return 0;
    }
}
