#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

typedef struct { unsigned char Cij[(n_par - 1) * d_par]; double fitness;} individual;
double uniform(void); /* function from stackoverflow, uses urandom file which stores entropy of the machine */
double ran1(long *idum);
void initializePopulation(individual** population, int popsize, fitness_function f);
void mutation(individual *ind, double probmutation);
void OnePointCrossover(individual* indivp1, individual* indivp2,
                       individual* indivf1, individual* indivf2);
individual* TournamentSelection(individual **population, int popsize);
void findFittest(individual* fittest, individual* population, int popsize);
int geneticAlgorithm(int popsize, double probmutation, fitness_function f, write_function writesol);
