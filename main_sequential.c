#include "myProtoSeq.h"


/*The main function is the entry point of the program. It reads input data from the "input.txt" file, 
calculates proximity points for different t-values, writes the results to the "output_seq.txt" file, 
and measures the execution time of the sequential implementation.
Return: An integer representing the exit status of the program.*/
int main(int argc, char *argv[])
{
    //Start the clock to measure the execution time of the program
    clock_t startTime, endTime;
    startTime = clock();
    //Declare variables for storing points data, N, K, D, tCount, tValues array, and proximites array
    Point *allPoints = NULL;
    int N, K, tCount;
    double D;
    double *tValues; // actualTs
    int *proximites; // Array that store proximity criteria

    printf("\n");
    //Read input data from the "input.txt" file and store it in appropriate variables using the readFromFile function
    allPoints = readFromFile(&N, &K, &D, &tCount); /*Reading from the file*/
    //Allocate memory for the tValues array to store the calculated t-values using calcTValues function
    tValues = (double *)malloc(sizeof(double) * (tCount + 1));
    if (!tValues)
    {
        fprintf(stderr, "Cannot Allocate memory.\n");
        exit(1);
    }
    calcTValues(tCount, tValues); /*This function will build the Tcount Array*/

    //Allocate memory for the proximites array to store the proximity points for each t-value and point using calculateProximity function
    proximites = (int *)malloc((tCount+1) * CONSTRAINTS * sizeof(int));
    if (!proximites)
    {
        fprintf(stderr, "Failed to allocate proximites array.\n");
        exit(1);
    }
    //// Initialize the array with -1 values
    for (int i = 0; i < (tCount+1) * CONSTRAINTS; i++) 
    {
        proximites[i] = -1;
    }
    calculateProximity(allPoints, N, tValues, tCount, D, proximites, K); /*Finding if the points are proimity criteria*/
    //Write the results to the "output_seq.txt" file using the writeOutputFile function
    writeOutputFile("output_seq.txt", tValues, tCount, proximites, allPoints, N);
    //End the clock to measure the execution time and calculate the elapsed time
    endTime = clock();
    double res = ((double) endTime - startTime) / CLOCKS_PER_SEC;
    printf("Sequential - Execution time: %.4lf seconds\n", res);

    return 0;
}
