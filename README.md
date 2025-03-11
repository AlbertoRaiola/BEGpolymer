# BEGpolymer
## Monte Carlo simulations of the Blume - Emery - Griffith magnetic polymer  

This code allows to perform Monte Carlo simulations of a Blume-Emery-Griffith magnetic polymer.  
A Monte Carlo sweep involves a global moves (pivot), N/4 local moves and N Glauber spin flip moves (N is the degree of polymerization).   
Multiple Markov chain (Parallel tempering) allows to exchange configurations at nearest neighbors temperatures.  
In the following a short desription of how the code work and its outputs and a brief description of each file.

## What this code do

- Variables and hash tables set up;
- Relaxation phase: the system undergoes a certain number of MC sweeps in order to termalize;
- Production phase: statistics regarding energy, gyration radius, magnetization, number of vacancies and number of contacts is collected.

Many quantities (such as specific heat, magnetic susceptibility, contacts variance, ...) are computed during the simulations and printed in the .tab file. 
The same occurs for autocorrelation times (collected in the .dia file). Moreover, time series of energy, magnetization, contacts, gyration radius and vacancies are collected in the .csv files. 

The list of input quantities can be found at rows 139-210 of the isaw3d_h.c file.  
An example of input file is test.in.

Output files are:

- datiE.csv datiC.csv, datiX.csv, datiRG.csv, datiM.csv (time series of observables);
- t400d2_*.dat (a file for each MC chain, where configurations are printed);
- t400d7.tab (statistics, in particular the first column is the inverse temperature and the height column is the specific heat);
- t400d2.dia (energy, magnetization, gyration radius and vacancies with error, autocorrelation times and their errors);
- t400d2.tab (garbage).

## Brief description of each routine

- allo.c and array_alloc.c: routine neeeded to allocate dynamically vectors/matrices/tensors of different dimensions;
- hash_chain.c: definition and functions for the hash_table (query the table, add and remove vertices);
- marsagliazo.c: generation of good sequences of (pseudo)random numbers;
- ising_saw_3d_h.c: spin flips Glauber moves;
- local_chain_ete.c: perform local moves on the polymeric substrate;
- pivot_chain_ete.c: perform pivot global moves on the polymeric substrate (needed to preserve ergodicity, however at very low temperature it is difficult for them to be accepted);
- isaw3d_h.c: main file

## How to compile the code

If you want to compile the code download it and put it in a folder. Open your shell it that folder and type:

cc -O3 *.c -o code.out -lm

## How to test the code

It is a good idea to put the input data on a file, like the test.in file in this repository. Type:

./code.out < test.in >> log &

In this case the statistics for the N = 50, Delta = 6 and K/J = 3 case is obtained. 













