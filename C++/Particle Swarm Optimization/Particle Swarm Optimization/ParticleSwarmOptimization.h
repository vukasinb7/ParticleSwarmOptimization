#include <tuple>
using namespace std;
using namespace tbb;
// Generate next iteration coordinates
// Param: start row, end row, number of iteration
void SerialNewPointGenerator(int iStart, int iEnd, int iter);
void ParallelNewPointGenerator(int iStart, int iEnd, int iter);
// Calculate initial function value and sets global best value
// Param: start index, end index, array
void SerialFunctionFirstCalculation(int iStart, int iEnd, double* arr, double (*func)(double* arr, int startIndex));
void ParallelFunctionFirstCalculation(int iStart, int iEnd, double* arr, double (*func)(double* arr, int startIndex));

// Calculate function value
// Param: start index, end index, array
void SerialFunctionCalculation(int iStart, int iEnd, double* arr, double (*func)(double* arr, int startIndex));
void ParallelFunctionCalculation(int iStart, int iEnd, double* arr, double (*func)(double* arr, int startIndex));

// Finds global best value
// Param: start index, end index
void SerialMinimum(int iStart, int iEnd, double *arr);
void ParallelMinimum(int iStart, int iEnd, double*);

// Particle Swarm Optimization Algorithm
tuple<double, double*> PSOParallel(double (*func)(double* arr, int startIndex), int N);
tuple<double, double*> PSOSerial(double (*func)(double* arr, int startIndex), int N);


