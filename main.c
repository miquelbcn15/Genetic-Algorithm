#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "utils.h"
#include "genetics.h"
#include "RKF78.h"

int main(int argc, char* argv[]) {
    // long idum = 0;
    // individual individual1;
    // individual individual2;
    // int randomData;
    // if ( (randomData = open("/dev/urandom", O_RDONLY)) < 0) {
    //         ExitError("while opening urandom", 1);
    // }
    // ssize_t result;
    // if ( (result = read(randomData, individual1.Cij, (n_par - 1)*d_par)) < 0) {
    //     ExitError("while generating random number", 1);
    // }
    // for (int i = 0; i < n_par - 1; i++) {
    //     for (int j = 0; j < d_par; j++) {
    //         individual1.Cij[i * d_par + j] = individual1.Cij[i * d_par + j] >> 4;
    //     }
    // }
    // for (int i = 0; i < n_par - 1; i++) {
    //     for (int j = 0; j < d_par; j++) {
    //         fprintf(stdout, "C[%d][%d] : %2d\n", i, j, (((int)((255) * ran1(&idum)) >> 4)));
    //     }
    // }
    // for (int i = 0; i < 100; i++) {
    //     int random = uniform() * ((n_par - 1) * d_par - 1);
    //     int posCoef = random / 4;
    //     int posBit  = random % 4;
    //     individual1.Cij[posCoef] = individual1.Cij[posCoef] ^ (1U << posBit);
    //     fprintf(stdout, "modification: %2d \n", individual1.Cij[posCoef]);
    // }

    int popsize;
    double probmut;

    if (argc<3
        || sscanf(argv[1], "%d", &popsize)!=1
        || sscanf(argv[2], "%lf", &probmut)!=1
       ) {
        fprintf(stderr,"%s popsize probmut \n", argv[0]);
        return -1;
    }
    geneticAlgorithm(popsize, probmut, &Curative_Fitness);

    return 0;
}
