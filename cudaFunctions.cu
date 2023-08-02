#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "myProto.h"



/*performing a smooth oscillatory motion between two points (x1 and x2) over time (t). 
It calculates the position x at a given time t using the sine function, which results in a periodic 
back-and-forth movement between the two points. */
__device__ double calculateX(double x1, double x2, double t)
{
/*
*param x1 (double): The initial position of the first point.
*param x2 (double): The initial position of the second point.
*param t (double): The time value for which the x value is to be calculated.
**/
    return ((x2 - x1) / 2) * sin(t * M_PI / 2) + ((x2 + x1) / 2);
}


//performing a linear transformation of the 'x' value to calculate the corresponding 'y' value based on the equation y = ax + b. 
__device__ double calculateY(double a, double b, double x)
{
/*
*param a (double): The coefficient 'a' in the linear equation y = ax + b.
*param b (double): The constant term 'b' in the linear equation y = ax + b.
*param x (double): The 'x' value for which the 'y' value is to be calculated.
**/
    return a * x + b;
}


/*
computing the Euclidean distance between two points (p1 and p2) in a 2D space at a given time t. 
It combines calculations from the calculateX and calculateY functions to determine the positions of the points 
and then calculates the distance between them using the Euclidean distance formula
*/
__device__ double calculateDistance(const Point *p1, const Point *p2, double t)
{
/*
*param p1 (const Point*): A pointer to the first Point structure representing the initial state of the first point.
*param p2 (const Point*): A pointer to the second Point structure representing the initial state of the second point.
*param t (double): The time value for which the distance between the points is to be calculated.
**/
    double x1 = calculateX(p1->x1, p1->x2, t);
    double y1 = calculateY(p1->a, p1->b, x1);

    double x2 = calculateX(p2->x1, p2->x2, t);
    double y2 = calculateY(p2->a, p2->b, x2);

    //Calculate the Euclidean distance between the two points (p1 and p2) at time t.
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}



//This device function checks whether the proximity criteria between two points p1 and p2 are met at a given time t and distance threshold D.
//Return: A boolean value indicating whether the proximity criteria are met.
__device__ bool isProximityCriteriaMet(const Point *p1, const Point *p2, double t, double D)
{
    /**
    *param p1: A pointer to the first Point structure representing one of the points.
    *param p2: A pointer to the second Point structure representing the other point.
    *param t: A double representing the time value.
    *param D: A double representing the distance threshold.
    **/

    //Get the id of the first point (p1) to avoid comparing the point with itself
    int currentPointIndex = p1->id;
  
    //Check if the points are the same (id of p1 and p2 are equal). 
    //If so, return false since a point is not considered in proximity to itself
    if (currentPointIndex == p2->id)
        return false;
    //Calculate the distance between the two points using the calculateDistance function
    double distance = calculateDistance(p1, p2, t);
    //Check if the calculated distance is less than or equal to the threshold D. If it is, return true, 
    //indicating that the proximity criteria are met. Otherwise, return false
    return distance <= D;
}



/*This device function updates the results array with the proximityPointId value for a specific time point index idx and constraint index
To avoid the race condition, you can use atomic operations to perform a compare-and-swap (CAS) operation. 
The CAS operation will only update the element if its value is still equal to -1, 
ensuring that only one thread successfully updates the result for each time point and constraint index.*/
__device__ void updateResults(int idx, int *results, int proximityPointId)
{
    /**
    *param idx: An integer representing the time point index.
    *param results: A pointer to the results array.
    *param proximityPointId: An integer representing the ID of the point that satisfies the proximity 
    **/

    //Loop through the constraints (indexed by j)
    for (int j = 0; j < CONSTRAINTS; j++)
    {
        // Calculate the index of the element to update in the results array based on the current idx and the constraint index j.
        int targetIndex = idx * CONSTRAINTS + j;
        // Use an atomic compare-and-swap (CAS) operation to update the value at targetIndex in the results array.
        // The CAS operation ensures that only one thread successfully updates the result for each time point and constraint index.
        int oldValue = -1;
        int newValue = proximityPointId;
        if (atomicCAS(&results[targetIndex], oldValue, newValue) == oldValue)
        {
            // The CAS operation succeeded, meaning this thread has updated the result successfully.
            // We can break out of the loop to avoid unnecessary checks.
            break;
        }
    }
}


/*This device function calculates the number of points in proximity to the point at currentPointIndex at a given time t and distance threshold D.
Return: An integer representing the number of points in proximity*/
__device__ int countProximityPoints(const Point *points, int N, int currentPointIndex, double t, double D)
{
/*
*param points (const Point*): A pointer to the array of Point structures representing the initial state of the points.
*param N (int): The total number of points in the array.
*param currentPointIndex (int): The index of the point for which proximity points are to be counted.
*param t (double): The time value for which the proximity criteria are to be checked.
*param D (double): The threshold distance used as the proximity criteria.
**/
    //Initialize a counter variable count to 0
    int count = 0;
    //Loop through all the points (j) to check their proximity to the point at currentPointIndex
    for (int j = 0; j < N; j++)
    {
        //Check if the current point (currentPointIndex) is different from the point being checked (j) to avoid comparing the point with itself.
        //Call the isProximityCriteriaMet function to determine if the proximity criteria are met between the points at 
        //currentPointIndex and j at time t and threshold distance D.
        if (currentPointIndex != j && isProximityCriteriaMet(&points[currentPointIndex], &points[j], t, D))
        {
            //If the proximity criteria are met, increment the count variable.
            count++;
        }
    }
    return count;
}


/*
checking the proximity of points in the points array at different time points represented by the tValues array. 
For each time point, the function iterates over all points and determines the number of points in proximity to each point based on the 
specified threshold distance D. If the number of proximity points (count) for a specific point is greater than or equal to the minimum 
required proximity points K, the point is considered a proximity point, and the results array is updated accordingly using the updateResults function.
*/
__global__ void checkProximityAtT(Point *points, double *tValues, const int tCount, const int N, const int K, const double D, int *results)
{
/*
*param points (Point*): A pointer to the array of Point structures representing the initial state of the points.
*param tValues (double*): A pointer to the array containing different time points.
*param tCount (const int): The total number of time points in the tValues array.
*param N (const int): The total number of points in the points array.
*param K (const int): The minimum number of proximity points required for a point to be considered a proximity point.
*param D (const double): The threshold distance used as the proximity criteria.
*param results (int*): A pointer to the results array.
**/
    //Calculate the index of the current thread
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    //Check if the calculated index is beyond the valid range of time points (tCount). -> do nothing
    if (idx >= tCount)
        return;

    //Retrieve the time value (t) corresponding to the current index from the tValues array.
    double t = tValues[idx];
    // Counter to track the number of points found for this t.
    int pointsFound = 0;

    for (int i = 0; i < N; i++)
    {
        //Call the countProximityPoints function to determine the number of points in proximity to the point at index i 
        //at time t and within the threshold distance D.
        int count = countProximityPoints(points, N, i, t, D);
        //Check if the count of proximity points (count) is greater than or equal to the minimum required proximity points K.
        if (count >= K)
        {
            //Retrieve the ID of the point at index i from the points array.
            int proximityPointId = points[i].id;
            // update the results array with the proximityPointId for the current time point index idx if the proximity criteria are met for the point at index i.
            updateResults(idx, results, proximityPointId);
            // Increment the counter for points found.
            pointsFound++; 
            // Break the loop if three points are found.
            if (pointsFound == 3)
                break;
        }
    }
}



//This function checks for any CUDA errors during GPU memory allocation and computation.
void checkCudaError(cudaError_t error, Point *dPoints, double *dTValues, int *dResults)
{
/**
*param error: A cudaError_t variable representing the CUDA error status.
*param dPoints: A pointer to the GPU memory for points data.
*param dTValues: A pointer to the GPU memory for tValues data.
*param dResults: A pointer to the GPU memory for results data.
**/
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA Error: %s\n", cudaGetErrorString(error));
        cudaFree(dPoints);
        cudaFree(dTValues);
        cudaFree(dResults);
        exit(EXIT_FAILURE);
    }
}

//This function performs the proximity computation on the GPU using CUDA
void performGPUComputation(int *N, int *K, double *D, int *tCountSize, double *myTValues, Point *points, int *results)
{
    //cudaSetDevice(0);
    
    cudaError_t error = cudaSuccess;
    //Calculate the number of threads per block for the GPU kernel computation, ensuring that the number of threads does not exceed BLOCK_SIZE.
    int threadPerBlock = min(BLOCK_SIZE, *tCountSize);
    //Calculate the number of blocks per grid for the GPU kernel computation based on the tCountSize.
    int blocksPerGrid = (*tCountSize + threadPerBlock - 1) / threadPerBlock;

    //Declare device pointers (dPoints, dTValues, dResults) to hold the corresponding data on the GPU.
    Point *dPoints = NULL;
    double *dTValues = NULL;
    int *dResults = NULL;

    //Allocate memory on the device for dPoints, dTValues, and dResults.
    error = cudaMalloc((void **)&dPoints, (*N) * sizeof(Point));
    checkCudaError(error, dPoints, dTValues, dResults);
    error = cudaMalloc((void **)&dTValues, (*tCountSize) * sizeof(double));
    checkCudaError(error, dPoints, dTValues, dResults);
    error = cudaMalloc((void **)&dResults, CONSTRAINTS * (*tCountSize) * sizeof(int));
    checkCudaError(error, dPoints, dTValues, dResults);

    //Copy the points data from the host to the device (dPoints) to enable GPU computations.
    error = cudaMemcpy(dPoints, points, (*N) * sizeof(Point), cudaMemcpyHostToDevice);
    checkCudaError(error, dPoints, dTValues, dResults);
    //Copy the tValues data from the host to the device (dTValues) to enable GPU computations.
    error = cudaMemcpy(dTValues, myTValues, (*tCountSize) * sizeof(double), cudaMemcpyHostToDevice);
    checkCudaError(error, dPoints, dTValues, dResults);
    //Copy the results data from the host to the device (dResults) to enable GPU computations.
    error = cudaMemcpy(dResults, results, CONSTRAINTS * (*tCountSize) * sizeof(int), cudaMemcpyHostToDevice);
    checkCudaError(error, dPoints, dTValues, dResults);

    // Launch the GPU kernel function checkProximityAtT with the specified grid and block configurations to perform 
    // proximity computation on the GPU.
    checkProximityAtT<<<blocksPerGrid, threadPerBlock>>>(dPoints, dTValues, *tCountSize, *N, *K, *D, dResults);

    //Synchronize the device, ensuring all GPU computations are completed before proceeding.
    error = cudaDeviceSynchronize();
    checkCudaError(error, dPoints, dTValues, dResults);

    //Copy the results of the GPU computation from the device (GPU) memory to the host (CPU) memory.
    error = cudaMemcpy(results, dResults, CONSTRAINTS * (*tCountSize) * sizeof(int), cudaMemcpyDeviceToHost);
    checkCudaError(error, dPoints, dTValues, dResults);

    //Release the device memory that was allocated for dPoints, dTValues, and dResults. This step is crucial to 
    //avoid memory leaks and efficiently manage GPU resources.
    cudaFree(dPoints);
    cudaFree(dTValues);
    cudaFree(dResults);
}

