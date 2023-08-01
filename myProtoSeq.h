#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#ifndef MYPROTO_H
#define MYPROTO_H

// Define constant values
#define CONSTRAINTS 3    // The number of constraints to satisfy

typedef struct {
    int id;
    double x1, x2, a, b;
} Point;


/**
* @param N: Pointer to store the value of N (number of points).
* @param K: Pointer to store the value of K (proximity criteria).
* @param D: Pointer to store the value of D (distance threshold).
* @param tCount: Pointer to store the total number of t-values.
* Return: A pointer to the array of Point structures.
**/
Point *readFromFile(int *N, int *K, double *D, int *tCount);


/**
 * @param tCount: The total number of t-values.
 * @paramtValues: Pointer to an array that will store the t-values.
 **/
void calcTValues(int tCount, double *tValues);


/**
* @param allPoints: Pointer to an array of all points (Point structures).
* @param N: The total number of points.
* @param tValues: Array of t-values.
* @param tCount: The total number of t-values.
* @param distance: The distance threshold for the proximity criteria.
* @param proximites: Pointer to an array that stores the proximity points.
* @param K: The number of points required to satisfy the proximity criteria.
**/
int calculateProximity(Point *allPoints, int N, double *tValues, int tCount, double distance, int *proximites, int K);

/**
* @param filename: The name of the output file to be written.
* @param tValues: Array of t-values.
* @param tCount: The total number of t-values.
* @param proximite: Array containing the proximity points.
* @param points: Array of Point structures.
* @param N: The total number of points.
**/
void writeOutputFile(const char *filename, double *tValues, int tCount, int *proximite, Point *points, int N);


#endif // MYPROTO_H




