#MPICC=vtcxx -vt:cxx mpic++
MPICC=mpic++
CC=g++
CCFLAGS= -O3 -Wall 
LDFLAGS=-lgsl -lgslcblas 


TARGET1=graphPart
TARGET2=isingTest
OBJS1= graph.o graphPartitioning.o sparse.o sparseBinary.o
OBJS2= graph.o distribGraph.o ising.o distribIsing.o sparse.o sparseBinary.o
HEADERS1=Graph.hpp GraphFactory.hpp Sparse.hpp SparseBinary.hpp GraphPartitioning.hpp Runlog.hpp RngWrapper.hpp
HEADERS1=Graph.hpp GraphFactory.hpp Sparse.hpp SparseBinary.hpp Runlog.hpp DistribGraph.hpp Ising.hpp RngWrapper.hpp


all: isingTest partTest

partTest: partTest.cc $(OBJS1) $(HEADERS1)
	$(CC) $(CCFLAGS) $(LDFLAGS) $(@).cc $(OBJS1) -o $(TARGET1)

isingTest: isingTest.cc $(OBJS2) $(HEADERS2)
	$(MPICC) $(CCFLAGS) $(LDFLAGS) $(@).cc $(OBJS2) -o $(TARGET2)

graph.o: graph.cc GraphFactory.hpp Sparse.hpp
	$(CC) -c $(CCFLAGS) graph.cc

graphPartitioning.o: graphPartitioning.cc GraphPartitioning.hpp
	$(CC) -c $(CCFLAGS) graphPartitioning.cc

ising.o: ising.cc Ising.hpp
	$(CC) -c $(CCFLAGS) ising.cc

distribGraph.o: distribGraph.cc DistribGraph.hpp Graph.hpp 
	$(MPICC) -c $(CCFLAGS) distribGraph.cc

distribIsing.o: distribIsing.cc DistribIsing.hpp 
	$(MPICC) -c $(CCFLAGS) distribIsing.cc

#checkEigen.o: checkEigen.cc CheckEigen.hpp
#	$(CC) -c $(CCFLAGS) checkEigen.cc

sparse.o: sparse.cc Sparse.hpp
	$(CC) -c $(CCFLAGS) sparse.cc

sparseBinary.o: sparseBinary.cc SparseBinary.hpp
	$(CC) -c $(CCFLAGS) sparseBinary.cc



dense: denseMesh.cc
	g++ denseMesh.cc -o dense

torus: torusMesh.cc
	g++ torusMesh.cc -o torus 

grid2d: rectMesh.cc
	g++ rectMesh.cc -o grid2d

randBivar: randBivar.cc
	g++ $(LDFLAGS) randBivar.cc -o randBivar

checkerMesh: checkerMesh.cc
	g++ checkerMesh.cc -o checkerMesh

testA2A: testAll2All.cc 
	$(MPICC) -c $(CCFLAGS) testAll2All.cc -o testA2A

test16:
	mpirun -np 16 ./isingTest paramsIsing.txt

test8:
	mpirun -np 8 ./isingTest paramsIsing.txt


clean:
	rm -f isingTest graphPart *.o