# Graph Partitioning for Computing on a (multi-core) Distributed Graph using MPI

Some code I wrote to perform Spectral Partitioning of Simple Graphs, and then, with the computed partition, perform distributed algorithms (currently an Ising Model) on the graph to test the timings.

Coded in C++ and MPI. The Graph parts are each stored locally on their respective nodes in a 'Distributed Graph' class which handles the communication across the Graph during computation.

###  Graph Partitioning

In this case done using Spectral Analysis (not the best, but very elegent theory).

Spectral analysis of a graph involves finding the Second (Fiedler) Eigenvector of the Lagrangian Matrix of the Graph.   

So I used the Lanczos Algorithm to transform the Langrange Matrix into a tridiagonal matrix, used Gerschgorin’s Theorem to find global bounds of Eigenvectors, a Sturm Sequence to find bounds on 2nd Eigenvector, then the Bisection Method to search the Characteristic Equation, the Inverse Power Method using LU decomposition to find the 2nd Eigenvector, and finally using the previously computed Lanczos Algorithms Transformation matrix Q, get the 2nd Eigenvector of the original matrix. 
This computation is so involved as it is for very large matrices, and targeting one Eigenvector using a Krylov Subspace method is the most practical route, maybe. 

### Distributed Computing 

Basically a Distributed Class which is a type of Class Sparse-Adjacency-Matrix-Graph, which is the baby of a Class Graph and Class Sparse-Binary-Matrix (Graph is a virtual class which needs a storage Class to become concrete). I built it all from scratch, even Sparse Matrix Type, for the exercise as it was for class.  
The algorithm I used to test the Distributed Graph was a Class Distributed-Ising,  which inherited from Class Ising, trying to keep a parallel between the serial and parallel verisons of my code. It uses Red-Green updating, so the graph needs to be a bipartite, very important restriction!  The Ising Model is calculated by MCMC/Gibbs sampling. 

MPI is built into the Distributed Classes but it is the Class Distributed-Graph that sets up the communicator that any algorithm would rely on.

Included are a number of toy Bipartite Graph examples (Torus-shaped Graph, Checkerboard Graph (with periodic boundary conditions) and Random Graph) to test the partitioning. Also a couple from Stanford's Large Network Dataset Collection.

There is a wrapper for the Random Number Generator so that another one can be swapped in easily, creating random numbers in a distributed way is handled by the algorithm.

