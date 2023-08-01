#include "myProtoSeq.h"
#include <string.h>


/*This function reads input data from the "input.txt" file and stores it in appropriate variables. 
It also reads the points data and returns an array of Point structures.*/
Point *readFromFile(int *N, int *K, double *D, int *tCount)
{
    Point *allPoints = NULL;
    //Open the "input.txt" file in read mode
    FILE *input_file = fopen("input.txt", "r");
    //Check if the file was successfully opened. If not, print an error message and exit the program
    if (input_file == NULL)
    {
        fprintf(stderr, "Failed to open input file.\n");
        exit(1);
    }
    //Read the first line of the file to get the values of N, K, D, and tCount
    if (fscanf(input_file, "%d %d %lf %d\n", N, K, D, tCount) != 4)
    {
        fprintf(stderr, "Failed reading from file.\n");
        fclose(input_file);
        exit(1);
    }
    //Allocate memory for the allPoints array to store the points data
    allPoints = (Point *)malloc((*N) * sizeof(Point));
    //Check if the memory allocation was successful. If not, print an error message, close the file, and exit the program
    if (!allPoints)
    {
        fprintf(stderr, "Cannot Allocate memory.\n");
        fclose(input_file);
        exit(1);
    }
    //Read N lines of (X1, X2, a, b) data for each point from the file and store them in the allPoints array
    for (int i = 0; i < *N; i++)
    {
        if (fscanf(input_file, "%d %lf %lf %lf %lf", &allPoints[i].id, &allPoints[i].x1, &allPoints[i].x2, &allPoints[i].a, &allPoints[i].b) != 5)
        {
            fprintf(stderr, "Failed reading from file.\n");
            fclose(input_file);
            free(allPoints);
            exit(1);
        }
    }
    //Close the input file and return the pointer to the array of Point structures
    fclose(input_file);
    return allPoints;
}

//This function calculates an array of t-values using a given tCount and stores them in the tValues array
/**
* @param tCount: The total number of t-values.
* @param tValues: Pointer to an array that will store the t-values.
**/
void calcTValues(int tCount, double *tValues)
{
    //Loop over each t-value and calculate it using the given formula
    for (int i = 0; i <= tCount; ++i)
        tValues[i] = (2.0 * i / tCount) - 1;
}


//This function updates the proximites array with a new point ID at the appropriate index
/**
* @param startingIndex: The starting index where the point ID should be updated in the proximites array.
* @param proximites: Pointer to an array that stores the proximity points.
* @param pointId: The ID of the point to be added to the proximites array.
**/
void updateProximitePoints(int startingIndex, int *proximites, int pointId) 
{
    //Loop over CONSTRAINTS number of times and check if there is an empty slot in the proximites array
    for (int i = 0; i < CONSTRAINTS; i++)
    {
        int index = startingIndex * CONSTRAINTS + i;
        //If an empty slot is found (indicated by the value -1), update it with the pointId and return from the function
        if (proximites[index] == -1)
        {
            proximites[index] = pointId; /*Put the point*/
            return;
        }
    }
}


//This function calculates the distance between two points at a given t value
/**
* @param p1: The first point (Point structure).
* @param p2: The second point (Point structure).
* @param t: The t-value at which the distance should be calculated.
* Return: The calculated distance between the two points.
**/
double calcDistance(const Point p1, const Point p2, double t)
{
    //Calculate the x and y coordinates for both points using the given formula
    double x1 = ((p1.x2 - p1.x1) / 2) * sin(t * M_PI / 2) + ((p1.x2 + p1.x1) / 2);
    double y1 = p1.a * x1 + p1.b;

    double x2 = ((p2.x2 - p2.x1) / 2) * sin(t * M_PI / 2) + (p2.x2 + p2.x1) / 2;
    double y2 = p2.a * x2 + p2.b;
    //Calculate the Euclidean distance between the two points and return the result
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}


//This function calculates the proximity points for each t-value based on the given distance and updates the proximites array accordingly.
int calculateProximity(Point *allPoints, int N, double *tValues, int tCount, double distance, int *proximites, int K)
{ 
    //nitialize a variable counter to keep track of the number of points that satisfy the proximity criteria at each t-value.
    int counter;
    //Loop over each t-value
    for (int i = 0; i <= tCount; i++) /*running on tCounts*/
    {
        //Loop over each point
        for (int j = 0; j < N; j++) /*Running on all the points*/
        {
            //Initialize counter to 0 for each t-value and point
            counter = 0;
            //Loop over each other point (except the current point) and check if the distance between them satisfies the proximity criteria.
            for (int k = 0; k < N; k++) /*check with other points but not the same point!*/
            {
                //If the distance satisfies the proximity criteria, increment the counter.
                if (j != k && calcDistance(allPoints[j], allPoints[k], tValues[i]) < distance)
                {
                    counter++;
                    //If the required number of points K is found to satisfy the proximity criteria, update the proximites array.
                    if (counter == K)
                    {
                        int pointId = allPoints[j].id;                 /*This point is proximity*/
                        updateProximitePoints(i, proximites, pointId); /*0,proximites,point that have 3 points*/
                        break;
                    }
                }
            }
        }
    }
}



//This function writes the proximity points that satisfy the criteria to the "output_seq.txt" file.
void writeOutputFile(const char *filename, double *tValues, int tCount, int *proximite, Point *points, int N)
{
    //Open the "output_seq.txt" file in write mode
    FILE *output_file = fopen("output_seq.txt", "w");
    //Check if the file was successfully opened. If not, print an error message and exit the program.
    if (!output_file)
    {
        fprintf(stderr, "Failed to open output file.\n");
        exit(1);
    }
    //Initialize a variable proximityFound to keep track of whether any set of points satisfying the proximity criteria is found or not.
    int proximityFound = 0;
    //Loop over each t-value and check if any set of points satisfies the proximity criteria.
    for (int i = 0; i <= tCount; i++)
    {
        //Initialize a variable counter to count the number of points that satisfy the proximity criteria for each t-value
        int counter = 0;
        //Calculate the starting index for the proximite array based on the current t-value
        int startIndex = CONSTRAINTS * i;
        //Loop over the CONSTRAINTS number of points for the current t-value and check if they satisfy the proximity criteria.
        for (int j = startIndex; j < startIndex + CONSTRAINTS; j++)
        {
            //If a point satisfies the proximity criteria (indicated by a non-negative value in proximite), increment the counter.
            if (proximite[j] != -1)
                counter++;
        }

        //If the required number of points CONSTRAINTS is found to satisfy the proximity criteria, 
        //update proximityFound and write the results to the output file
        if (counter == CONSTRAINTS)
        {
            proximityFound = 1;

            fprintf(output_file, "Points ");
            for (int j = startIndex; j < CONSTRAINTS + startIndex; j++)
            {
                fprintf(output_file, "pointID%d", proximite[j]);
                if (j < CONSTRAINTS + startIndex - 1)
                    fprintf(output_file, ", ");
            }
            fprintf(output_file, " satisfy Proximity Criteria at t = %.2f\n", tValues[i]);
        }
    }
    //If no points satisfy the proximity criteria, write a message to the output file
    if (!proximityFound)
        fprintf(output_file, "There were no 3 points found for any t.\n");
    //Close the output file
    fclose(output_file); 
}
