#pragma once
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// Define constant values
#define BLOCK_SIZE 256   // The block size used for GPU computations
#define CONSTRAINTS 3    // The number of constraints to satisfy

# ifndef POINT_H
# define POINT_H

// Define the Point structure
typedef struct Point {
    int id;             // Point ID
    double x1, x2, a, b; // Point attributes
} Point;

#endif


// Function prototypes

/**
* @param rank (int): The rank of the current process.
* @param size (int): size: The total number of processes.
* @param tCount (int): The total number of elements in the array to be distributed.
* @param sendcounts (int*): An array of size size, to be filled with the number of elements each process will receive.
* @param displs (int*): An array of size size, to be filled with the displacement (starting index) of each process in the global array.
**/
void calculateSendcountsAndDispls(int rank, int size, int tCount, int *sendcounts, int *displs);


/* *
* @param myTValuesSize (int): The size of the array to be allocated.
* Return: A pointer to the allocated memory block (the array).
**/
int *allocateResults(int myTValuesSize);


/*
* @param N (int*): A pointer to the number of points in the points array.
* @param K (int*): A pointer to the proximity constraint value K.
* @param D (double*): A pointer to the distance threshold D.
* @param tCountSize (int*): A pointer to the size of the myTValues array.
* @param myTValues (double*): A pointer to the array containing t values for the current process.
* @param points (Point*): A pointer to the array of Point structures representing all points.
* @param results (int*): A pointer to the results array that will store the proximity point IDs.
**/
void performGPUComputation(int *N, int *K, double *D, int *tCountSize, double *myTValues, Point *points, int *results);


/**
* @param filename (const char*): The name of the output file where the results will be written.
* @param tValues (double*): An array of t values used in the computation.
* @param tCount (int): The total number of t values in the tValues array.
* @param results (int*): An array containing the computed results from the computation.
* @param points (Point*): An array of Point structures representing data points (not directly used in this function).
* @param N (int): The total number of data points (not directly used in this function).
**/
void writeResultsToFile(const char* filename, double* tValues, int tCount, int* results, Point* points, int N);



/**
* @param rank (int): The rank or identifier of the current process within the group of processes.
* @param size (int): The total number of processes in the MPI communicator.
* @param N (int): The total number of data points (not directly used in this function).
* @param tCount (int): The total number of time values.
* @param tCountSize (int): The number of time values for the current process.
* @param results (int*): A pointer to the array containing the computed results for the current process's time values.
* @param globalResults (int*): A pointer to the array where the gathered results will be stored on the root process.
**/
void collectResults(int rank, int size, int N, int tCount, int tCountSize, int *results, int *global_results);



/**
* @param tCount (int): The number of elements in the tValues array to be calculated.
* @param tValues (double**): Pointer to a pointer that will store the array of tValues.
**/
void calculateTValues(int tCount, double **tValues);



/**
* @param filename (const char*): The name of the input file containing the data.
* @param N (int*): A pointer to an integer variable to store the number of data points read from the file.
* @param K (int*): A pointer to an integer variable to store a minimum number of points read from the file.
* @param D (double*): A pointer to a double variable to store another parameter read from the file.
* @param tCount (int*): A pointer to an integer variable to store the number of time values read from the file.
* @param points (Point**): A pointer to a pointer to the array of data points (structure representing a point in space).
**/
void readInputFromFile(const char *filename, int *N, int *K, double *D, int *tCount, Point **points);




