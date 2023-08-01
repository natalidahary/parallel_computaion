#include "myProto.h"
#include <mpi.h>
#include <string.h>

/*This function calculates the sendcounts and displacements for distributing an array of elements among different processes
 in a parallel environment using MPI (Message Passing Interface). The goal is to evenly distribute the elements across processes, 
 considering any remainders for a fair distribution.
 */
void calculateSendcountsAndDispls(int rank, int size, int tCount, int *sendcounts, int *displs)
{
    // Calculate the average number of elements that each process will receive (assuming an even distribution).
    int tCountSize = tCount / size;
    // Calculate the remaining number of elements that need to be distributed after an even distribution is performed.
    int remainTValues = tCount % size;

    //Set the base sendcount and displacement for the current process:
    int baseSendcount = tCountSize;
    int baseDisplacement = rank * tCountSize;

    //If the current process rank is less than remainTValues, it means it should receive an additional element, 
    //so increment the base sendcount and adjust the base displacement accordingly:
    if (rank < remainTValues) {
        baseSendcount++;
        baseDisplacement += rank;
    } else {
        // Processes with rank greater than or equal to remainTValues do not receive an additional element.
        baseDisplacement += remainTValues;
    }

    //Loop over all processes and set their sendcounts and displacements:
    for (int i = 0; i < size; i++) {
        // Set the sendcount for the current process.
        sendcounts[i] = baseSendcount;
        // Set the displacement for the current process.
        displs[i] = baseDisplacement;

        // Move the baseDisplacement to the next position for the next iteration.
        baseDisplacement += baseSendcount;

        // If the next process will receive an additional element, increment baseSendcount.
        if (i + 1 < remainTValues) {
            baseSendcount++;
        }
    }
}

/* This function dynamically allocates memory for an array that will store the results of a computation. 
It initializes all elements of the array to the value -1.
*/
int *allocateResults(int myTValuesSize)
{
    //Dynamically allocate memory for the results array using malloc:
    int *results = (int *)malloc(CONSTRAINTS * myTValuesSize * sizeof(int));
    //Use memset to set all elements of the results array to the value -1:
    memset(results, -1, CONSTRAINTS * myTValuesSize * sizeof(int));
    //Return the pointer to the allocated memory
    return results;
}


// Helper function to read values from the input file - N   K   D   TCount
static void readValues(FILE *file, int *value1, int *value2, double *value3, int *value4) {
    if (fscanf(file, "%d %d %lf %d\n", value1, value2, value3, value4) != 4) {
        fprintf(stderr, "Failed to read required values from input file.\n");
        fclose(file);
        MPI_Finalize();
        exit(1);
    }
}


// Helper function to read point data from the input file and populate the points array - id   x1    x2    a    b
static void readPointData(FILE *file, int N, Point **points) {
    for (int i = 0; i < N; ++i) {
        if (fscanf(file, "%d %lf %lf %lf %lf\n", &((*points)[i].id), &((*points)[i].x1), &((*points)[i].x2), &((*points)[i].a), &((*points)[i].b)) != 5) {
            fprintf(stderr, "Failed to read point data from input file.\n");
            fclose(file);
            MPI_Finalize();
            exit(1);
        }
    }
}


//reading input data from a specified file in a parallel computing environment using MPI
void readInputFromFile(const char *filename, int *N, int *K, double *D, int *tCount, Point **points) {
    // Open the input file in read mode
    FILE *file = fopen(filename, "r");
    //Check if the file was successfully opened. If not, print an error message and exit the program
    if (!file) {
        fprintf(stderr, "Failed to open input file.\n");
        MPI_Finalize();
        exit(1);
    }

    //Read the values N, K, D, and tCount from the file using the readValues function
    readValues(file, N, K, D, tCount);

    //Allocate memory for the points array using malloc, based on the value of N
    *points = (Point *)malloc(*N * sizeof(Point));
    //Check if the memory allocation was successful. If not, print an error message, close the file, and exit the program
    if (!(*points)) {
        fprintf(stderr, "Failed to allocate points.\n");
        fclose(file);
        MPI_Finalize();
        exit(1);
    }
    //Read the point data from the file using the readPointData function- id   x1    x2    a    b
    readPointData(file, *N, points);
    // Close the input file
    fclose(file); 
}


//This function calculates an array of t-values, which are used for some computation.
void calculateTValues(int tCount, double **tValues) {
    //Allocate memory for the t-values array using malloc
    *tValues = (double *)malloc((tCount + 1) * sizeof(double)); 

 //This OpenMP directives to parallelize the for loop, allowing multiple threads to efficiently compute the t values concurrently.
#pragma omp parallel for
    for (int i = 0; i <= tCount; ++i)
        //Calculate the t value for each index i and store it in the tValues array at the corresponding index i.
        (*tValues)[i] = (2.0 * i / tCount) - 1; 
}


/*
gathering the computed results from all processes into a global 2D array on the root process (rank 0). 
It uses MPI's MPI_Gatherv function, which allows efficient data gathering and placement of the results in the globalResults array. 
*/
void collectResults(int rank, int size, int N, int tCount, int tCountSize, int *results, int *globalResults) {
    //storing the number of elements that each process contributes to the gathered array.
    int *recvcounts = (int *)malloc(size * sizeof(int));
    //storing the displacements (offsets) for each process's data in the gathered array.
    int *displs = (int *)malloc(size * sizeof(int));

    //Calculate the recvcounts array to determine the number of elements each process will contribute to the gathered array.
    for (int i = 0; i < size; i++) {
        recvcounts[i] = CONSTRAINTS * ((i < tCount % size) ? (tCount / size + 1) : (tCount / size));
    }

    //Calculate the displs array to determine the displacement (starting index) for each process's data in the gathered array.
    displs[0] = 0;
    for (int i = 1; i < size; i++) {
        //The displacement is the sum of the previous process's displacement and the number of elements contributed
        //by the previous process (recvcounts[i - 1]).
        displs[i] = displs[i - 1] + recvcounts[i - 1];
    }

    //Gather the 2D array results from all processes into the globalResults array on the root process using the MPI function MPI_Gatherv.
    MPI_Gatherv(results, CONSTRAINTS * tCountSize, MPI_INT,
                globalResults, recvcounts, displs, MPI_INT,
                0, MPI_COMM_WORLD);
    //Free the memory allocated for the recvcounts and displs arrays.
    free(recvcounts); 
    free(displs);     
}



/*
 writing the computed results from a parallel computation to an output file named filename. 
 It checks if any set of points satisfying the proximity criteria is found for each t value and writes the results to the output file.
*/
void writeResultsToFile(const char *filename, double *tValues, int tCount, int *results, Point *points, int N) {
    //Open the output file specified by filename in write mode ("w").
    FILE *file = fopen(filename, "w");
    //Check if the file was successfully opened. If not, print an error message and exit the program
    if (!file) {
        fprintf(stderr, "Failed to open output file.\n");
        MPI_Finalize();
        exit(1);
    }

    //Initialize a variable proximityFound to keep track of whether any set of points satisfying the proximity criteria is found or not.
    int proximityFound = 0;

    //Loop over all t-values and check if any set of points satisfies the proximity criteria.
    for (int i = 0; i < tCount; i++) {
        int count = 0;
        //Initialize an array pointIDs to store the IDs of the points that satisfy the proximity criteria. 
        //The array is initialized with -1 values.
        int pointIDs[CONSTRAINTS] = {-1, -1, -1};

        for (int j = 0; j < CONSTRAINTS; j++) {
            //Get the ID of the point
            int pointID = results[i * CONSTRAINTS + j];
            // Check if the pointID is within the valid range of data points.
            if (pointID >= 0 && pointID < N) {
                //Store the valid pointID in the pointIDs array.
                pointIDs[count] = pointID;
                //Increment the count to track the number of valid points found for the current t value.
                count++;
            }
        }

        // Check if the number of valid points (count) is equal to the constraint value (CONSTRAINTS). 
        //If true, it means that the proximity criteria are satisfied.
        if (count == CONSTRAINTS) {
            //at least one set of points satisfying the proximity criteria has been found.
            proximityFound = 1;

            // Print the string "Points " to the output file.
            fprintf(file, "Points ");
            for (int j = 0; j < CONSTRAINTS; j++) {
                fprintf(file, "pointID%d", pointIDs[j]);
                if (j < 2)
                    fprintf(file, ", ");
            }
            fprintf(file, " satisfy Proximity Criteria at t = %.2f\n", tValues[i]);
        }
    }
    //Check if no set of points satisfying the proximity criteria was found for any t value.
    if (!proximityFound) {
        fprintf(file, "There were no 3 points found for any t.\n");
    }

    fclose(file); // Close the output file
}
