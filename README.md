# Genetic Algorithm's 

This code corresponds to the Genetic Algorithm assignment from the course of
Optimization taken during the Master in Modeling for Science and Engineering 
at the Autonomous University of Barcelona  2019-2020.
The program is based on the cancer modelling problem from [McC05, Example 5.3].
It uses the well-known Genetic Algorithm to find a minimum for a given fitness
solution proposed by Gompertz. If this is not possible, it tries to find a
paliative treatment.

## Code structure and how to run it

The program is structured into three main files:
    - utils.c: contain some given functions by Lluis Alseda and others focused 
    on writing the results.
    - genetics.c: contains the main functions corresponding to the Genetic
    Algorithm implementation.
    - main.c: includes the main function, it basically runs the code
    - The functions corresponding to the files RKF78.c and its header are
    developed and provided by Lluis Alseda.

### Instructions
The src folder contains a Makefile to compile the code. The main object file
receives two parameters by command line and another optional one. 
    - int probsize: size of the population of individuals.
    - double probmut: probability of mutation applied in the mutation method.
    - int typesol: (optional) type of solution strategy, by default it is 0.
    The options are the following ones:
        - 0: tries to find a curative solution, otherwise, it tries
        to find a paliative solution.
        - 1: it tries to find straightly a paliative solution.
        - 2: it tries to find botw solutions.
Depending on the type of solution used the code generates two main solution
files:
    - curative\_solution.sol 
    - paliative\_solution.sol
Both files contains the values of (t,N) and at the end the found concentrations by the Algorithm.
The program can be executed by running the script chou.sh.

```bash
bash chou.sh
