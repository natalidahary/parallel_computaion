NUM_PROC = 2
build:
	gcc ./main_sequential.c cFunctionsSeq.c -o mainSeq -lm
	mpicxx -fopenmp -c main_parallel.c -o main_parallel.o
	mpicxx -fopenmp -c cFunctions.c -o cFunctions.o
	nvcc -I./Common  -gencode arch=compute_61,code=sm_61 -c cudaFunctions.cu -o cudaFunctions.o
	mpicxx -fopenmp -o mpiCudaOpemMP  main_parallel.o cFunctions.o cudaFunctions.o  -L/usr/local/cuda/lib -L/usr/local/cuda/lib64 -lcudart
	

clean:
	rm -f mainSeq
	rm -f *.o ./mpiCudaOpemMP

run:	
	./mainSeq
	mpiexec -np ${NUM_PROC} ./mpiCudaOpemMP

runOn2:
	mpiexec -np ${NUM_PROC} -machinefile  mf  -map-by  node  ./mpiCudaOpemMP
