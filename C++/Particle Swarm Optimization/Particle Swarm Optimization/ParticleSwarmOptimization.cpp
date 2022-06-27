#include <iostream>
#include <vector>
#include <tuple>
#include <math.h>
#include <limits>
#include <algorithm>
#include <functional>
#include "tbb/tick_count.h"
#include "tbb/task_group.h"
#include "ParticleSwarmOptimization.h"

extern int DIMENSIONS;
extern int PARTICLE_NUM;
extern int N;

extern int CUT_OFF_PARTICLE_NUM;
extern const unsigned int DATA_SIZE;

extern double* X, * XStart; // particle coordinates
extern double* V, * VStart; // particle velocity

extern double* pBestX, * pBestXStart; //personal best function values for each particle
extern double* pBestValue, * pBestValueStart; // coordinates for personal best value for each particle

extern double gBestValue; //global best function value
extern double* gBestX, * gBestXStart; //coordinates for global best function value

extern double* rp; //random numbers
extern double* rg; //random numbers

#pragma region ParallelizePointGenerator
void SerialNewPointGenerator(int iStart, int iEnd, int iter)
{
	double w = 0.9 - iter / N * (0.9 - 0.4); //const of intertia
	double cp = 2.5 - iter / N * (2.5 - 0.5); //cognitive const
	double cg = 0.5 + iter / N * (2.5 - 0.5); //social const

	for (int i = iStart;i < iEnd;i++) {
		for (int j = 0;j < DIMENSIONS;j++) {
			V[i * DIMENSIONS + j] = V[i * DIMENSIONS + j] * w + ((X[i * DIMENSIONS + j] - 
				pBestX[i * DIMENSIONS + j]) * (-1)) * cp * rp[iter] + ((X[i * DIMENSIONS + j] - 
					gBestX[j]) * (-1)) * cg * rg[iter];//calculate new velocity
			X[i * DIMENSIONS + j] = X[i * DIMENSIONS + j] + V[i * DIMENSIONS + j]; //calculate new coordinates
		}
	}
}


void ParallelNewPointGenerator(int iStart, int iEnd, int iter) {
	task_group g;
	if ((iEnd - iStart) < CUT_OFF_PARTICLE_NUM) {
		SerialNewPointGenerator(iStart, iEnd, iter);
	}
	else {
		g.run([&] {ParallelNewPointGenerator(iStart, (iStart + iEnd) / 2, iter);});
		g.run([&] {ParallelNewPointGenerator((iStart + iEnd) / 2, iEnd, iter);});
		g.wait();
	}
}
#pragma endregion ParallelizePointGenerator

#pragma region ParallelizeFirstCalculation
void SerialFunctionFirstCalculation(int iStart, int iEnd, double* arr,
	double (*func)(double* arr, int startIndex)) {
	for (int j = iStart;j < iEnd;j++) {
		arr[j] = (func(X, j)); //calculate function value for new coordinates
		if (arr[j] < gBestValue) { //set best global function value
			gBestValue = arr[j];
			for (int k = 0;k < DIMENSIONS;k++) {
				gBestX[k] = pBestX[j * DIMENSIONS + k];
			}
		}
	}
}
void ParallelFunctionFirstCalculation(int iStart, int iEnd, double* arr, 
	double (*func)(double* arr, int startIndex)) {
	task_group g1;
	if ((iEnd - iStart) < CUT_OFF_PARTICLE_NUM) {
		SerialFunctionFirstCalculation(iStart, iEnd, arr,func);
	}
	else {
		g1.run([&] {ParallelFunctionFirstCalculation(iStart, (iStart + iEnd) / 2, arr,func);});
		g1.run([&] {ParallelFunctionFirstCalculation((iStart + iEnd) / 2, iEnd, arr,func);});
		g1.wait();

	}
}

#pragma endregion ParallelizeFirstCalculation


#pragma region ParallelizeOnlyCalculation
void SerialFunctionCalculation(int iStart, int iEnd, double* arr, double (*func)(double* arr, int startIndex)) {
	for (int j = iStart;j < iEnd;j++) {
		arr[j] = (func(X, j));
	}
}
void ParallelFunctionCalculation(int iStart, int iEnd, double* arr, double (*func)(double* arr, int startIndex)) {
	task_group g2;
	if ((iEnd - iStart) < CUT_OFF_PARTICLE_NUM) {
		SerialFunctionCalculation(iStart, iEnd, arr,func);
	}
	else {
		g2.run([&] {ParallelFunctionCalculation(iStart, (iStart + iEnd) / 2, arr,func);});
		g2.run([&] {ParallelFunctionCalculation((iStart + iEnd) / 2, iEnd, arr,func);});
		g2.wait();

	}
}

#pragma endregion ParallelizeOnlyCalculation

#pragma region ParallelizeFindingMinimum
void SerialMinimum(int iStart, int iEnd,double *arr) {
	for (int i = iStart;i < iEnd;i++) {
		if (arr[i] > pBestValue[i])
		{
			pBestValue[i] = arr[i];
			for (int j = 0; j < DIMENSIONS; j++)
			{
				pBestX[i * DIMENSIONS + j] = X[i * DIMENSIONS + j];
			}
			if (gBestValue > pBestValue[i]) {
				gBestValue = pBestValue[i];
				for (int k = 0;k < DIMENSIONS;k++) {
					gBestX[k] = pBestX[i * DIMENSIONS + k];
				}
			}
		}
		
	}
}
void ParallelMinimum(int iStart, int iEnd, double * arr) {
	task_group g3;
	if ((iEnd - iStart) < CUT_OFF_PARTICLE_NUM) {
		SerialMinimum(iStart, iEnd, arr);
	}
	else {
		g3.run([&] {ParallelMinimum(iStart, (iStart + iEnd) / 2, arr);});
		g3.run([&] {ParallelMinimum((iStart + iEnd) / 2, iEnd, arr);});
		g3.wait();

	}
}

#pragma endregion ParallelizeFindingMinimum



tuple<double, double*> PSOParallel(double (*func)(double* arr, int startIndex), int N) {

	ParallelFunctionFirstCalculation(0, PARTICLE_NUM, pBestValue, func);//initial calculation



	for (int i = 0;i < N;i++) {

		ParallelNewPointGenerator(0, PARTICLE_NUM, i);//generate new velocity and position


		double* newValue = new double[PARTICLE_NUM]; //new function values
		ParallelFunctionCalculation(0, PARTICLE_NUM, newValue, func);//calculate new values
		//set new personal and global best function value
		ParallelMinimum(0, PARTICLE_NUM,newValue);
		delete[] newValue;
	}
	return make_tuple(gBestValue, gBestX);
}

tuple<double, double*> PSOSerial(double (*func)(double* arr, int startIndex), int N) { 

	SerialFunctionFirstCalculation(0, PARTICLE_NUM, pBestValue, func);//initial calculation


	for (int i = 0;i < N;i++) {

		SerialNewPointGenerator(0, PARTICLE_NUM, i);//generate new velocity and position

		double* newValue = new double[PARTICLE_NUM]; //new function values
		SerialFunctionCalculation(0, PARTICLE_NUM, newValue, func);//calculate new values
		//set new personal and global best function value
		SerialMinimum(0, PARTICLE_NUM, newValue);
		delete[] newValue;
	}
	return make_tuple(gBestValue, gBestX);
}




