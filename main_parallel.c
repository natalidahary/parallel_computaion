#include <mpi.h>
#include <string.h>
#include "myProto.h"

/*The main function is the entry point of the parallel version of the program. It uses MPI (Message Passing Interface) to distribute 
the computation of proximity points across multiple processes. The input data is read from the "input.txt" file, 
and the results are written to the "output.txt" file. The function measures the execution time of the parallel implementation
Return: An integer representing the exit status of the program.
*/
int main(int argc, char *argv[])
{
    clock_t startTime, endTime;
    // Start the timer
    startTime = clock();
    // Initialize MPI environment
    MPI_Init(&argc, &argv);
    
    //Declare variables for storing points data, N, K, D, tCount, tValues array, proximites array, and MPI-related variables
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    char input_file[] = "input.txt";
    char output_file[] = "output.txt";
    int N, K, tCount;
    double D;
    Point *points = NULL;
    double *tValues = NULL;

    // Read input data on root process (rank 0)
    if (rank == 0)
    {
       readInputFromFile(input_file, &N, &K, &D, &tCount, &points);
       calculateTValues(tCount, &tValues);
    }
    //Broadcast the input data (N, K, D, tCount, tValues, and points) to all processes using MPI_Bcast
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&K, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&D, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Define the custom MPI datatype MPI_POINT to represent the Point structure.
    MPI_Datatype MPI_POINT;
    MPI_Type_contiguous(sizeof(Point), MPI_BYTE, &MPI_POINT);
    MPI_Type_commit(&MPI_POINT);

    // Allocate memory for tValues on non-root processes
    if (rank != 0)
    {
        tValues = (double *)malloc(tCount * sizeof(double));
        if (!tValues)
        {
            fprintf(stderr, "Failed to allocate tValues.\n");
            MPI_Finalize();
            exit(1);
        }
    }

    //Broadcast tValues from the root process to all processes using MPI_Bcast
    MPI_Bcast(tValues, tCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //Allocate memory for tValues on non-root processes
    if (rank != 0)
    {
        points = (Point *)malloc(N * sizeof(Point));
        if (!points)
        {
            fprintf(stderr, "Failed to allocate points.\n");
            MPI_Finalize();
            exit(1);
        }
    }
  

    //Broadcast points array from the root process to all processes using MPI_Bcast
    MPI_Bcast(points, N, MPI_POINT, 0, MPI_COMM_WORLD);

    // Calculate sendcounts and displs arrays for efficient distribution of time values among processes using the calculateSendcountsAndDispls function
    int *sendcounts = (int *)malloc(size * sizeof(int));
    int *displs = (int *)malloc(size * sizeof(int));
    calculateSendcountsAndDispls(rank, size, tCount, sendcounts, displs);

    // Determine the number of time values that the current process will receive
    int myTValuesSize = sendcounts[rank];
    // Allocate memory for the myTValues array to store the time values that the current process will work on
    double *myTValues = (double *)malloc(myTValuesSize * sizeof(double));

    //Scatter t values to all processes using MPI_Scatterv
    MPI_Scatterv(tValues, sendcounts, displs, MPI_DOUBLE, myTValues, myTValuesSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //Allocate memory for the results array using the allocateResults function
    int *results = allocateResults(myTValuesSize);

    //Compute results on the GPU using the performGPUComputation function
    performGPUComputation(&N, &K, &D, &myTValuesSize, myTValues, points, results);

    int *globalResults = NULL;
    if (rank == 0)
    {
        // Allocate memory for global results on root process
        globalResults = (int *)malloc(CONSTRAINTS * tCount * sizeof(int));
        // Initialize all elements to the default value (-1) to ensure proper data handling
        memset(globalResults, -1, CONSTRAINTS * tCount * sizeof(int));
    }

    //Gather results from all processes to the root process using the collectResults function
    collectResults(rank, size, N, tCount, myTValuesSize, results, globalResults);

    if (rank == 0)
    {
        // Write the final results to an output file named "output.txt" on the root process
        writeResultsToFile(output_file, tValues, tCount, globalResults, points, N);
           
       // Stop the timer before MPI_Finalize
       endTime = clock();
       double executionTime = ((double) endTime - startTime) / CLOCKS_PER_SEC;
       // Print execution time on rank 0
       printf("Parallel - Execution time: %.4lf seconds\n", executionTime); 
 
    }
    
    // Free memory on all processes
    free(globalResults);
    free(tValues);
    free(sendcounts);
    free(displs);
    free(myTValues);
    free(results);
    free(points);
    MPI_Type_free(&MPI_POINT);
 

    // Finalize MPI environment
    MPI_Finalize();

    return 0;
}


