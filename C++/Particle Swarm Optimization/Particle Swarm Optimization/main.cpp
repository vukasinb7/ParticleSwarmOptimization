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

using namespace std;
using namespace tbb;
int DIMENSIONS = 4;
int PARTICLE_NUM = 12000;
int N = 1000;

int CUT_OFF_PARTICLE_NUM = 800;
static const unsigned int DATA_SIZE = DIMENSIONS * PARTICLE_NUM;

double* X, * XStart; // koordinate cestica(nasumicno izabrane)
double* V, * VStart; // brzine cestica (nasumicno izabrane)

double* pBestX, * pBestXStart; //koordinate personalno najbolje vrednosti funkcije
double* pBestValue, * pBestValueStart; //personalno najbolja vrednost funkcije

double gBestValue = numeric_limits<double>::max(); //globalno najbolja vrednost funkcije
double* gBestX, * gBestXStart; //koordinate globalno najbolje funkcije.

double* rp;
double* rg;


double TestFunction(double x) {
	return pow(x, 6) - 4 * pow(x, 2) + 2 * pow(x, 4) + 24 * x;
}
double TestFunction2(double* x, int i) {
	return 0.01 * (pow((x[i * DIMENSIONS + 0] - 1), 2) + 2 * pow((x[i * DIMENSIONS + 1] - 1), 2)) *
		(pow((x[i * DIMENSIONS + 0] + 1), 2) + 2 * pow((x[i * DIMENSIONS + 1] + 1), 2) + 0.5) *
		(pow((x[i * DIMENSIONS + 0] + 2), 2) + 2 * pow((x[i * DIMENSIONS + 1] - 2), 2) + 0.7);
}
double TestFunction3(double* x, int i) {
	return 0.01 * (pow((x[i * DIMENSIONS + 2] - 1), 5) + 2 * pow((x[i * DIMENSIONS + 1] - 1), 4)) *
		(pow((x[i * DIMENSIONS + 0] + 1), 2) + 2 * pow((x[i * DIMENSIONS + 1] + 3), 2) + 0.5) *
		(pow((x[i * DIMENSIONS + 0] + 2), 2) + 2 * pow((x[i * DIMENSIONS + 0] - 2), 2) + 0.7) +
		0.01 * (pow((x[i * DIMENSIONS + 0] - 1), 2) + 2 * pow((x[i * DIMENSIONS + 3] - 1), 7)) *
		(pow((x[i * DIMENSIONS + 2] + 1), 1) + 2 * pow((x[i * DIMENSIONS + 1] + 3), 3) + 0.5) *
		(pow((x[i * DIMENSIONS + 1] + 2), 2) + 2 * pow((x[i * DIMENSIONS + 2] - 2), 3) + 0.7) + 2;
}

void Initialize() {
	XStart = new double[DATA_SIZE];
	VStart = new double[DATA_SIZE];
	pBestXStart = new double[DATA_SIZE];
	pBestValueStart = new double[PARTICLE_NUM];
	gBestXStart = new double[DIMENSIONS];

	X = new double[DATA_SIZE];
	V = new double[DATA_SIZE];
	pBestX = new double[DATA_SIZE];
	pBestValue = new double[PARTICLE_NUM];
	gBestX = new double[DIMENSIONS];
	rp = new double[N];
	rg = new double[N];
	for (int i = 0;i < PARTICLE_NUM;i++) {
		for (int j = 0; j < DIMENSIONS; j++) //postavljanje inicijalnih vrednosti
		{
			double x = (double)(rand()) / ((double)(RAND_MAX));
			double v = (double)(rand()) / ((double)(RAND_MAX));
			XStart[i * DIMENSIONS + j] = x;
			VStart[i * DIMENSIONS + j] = v;
			pBestXStart[i * DIMENSIONS + j] = x;
		}
	}
	for (int i = 0; i < N; i++)
	{
		rp[i] = (double)(rand()) / ((double)(RAND_MAX)); //random vrednosti izmedju 0 i 1
		rg[i] = (double)(rand()) / ((double)(RAND_MAX)); //random vrednosti izmedju 0 i 1
	}

}
void ResetParameters() {
	for (int i = 0;i < PARTICLE_NUM;i++) {
		for (int j = 0; j < DIMENSIONS; j++)
		{
			X[i * DIMENSIONS + j] = XStart[i * DIMENSIONS + j];
			V[i * DIMENSIONS + j] = VStart[i * DIMENSIONS + j];
			pBestX[i * DIMENSIONS + j] = pBestXStart[i * DIMENSIONS + j];
		}
		pBestValue[i] = pBestValueStart[i];


	}
	for (int i = 0; i < DIMENSIONS; i++)
	{
		gBestX[i] = gBestXStart[i];
	}
	gBestValue = numeric_limits<double>::max();
}
void Destruct() {
	delete[] X, XStart, V, VStart;
	delete[] pBestX,pBestXStart,pBestValue,pBestValueStart,gBestX,gBestXStart;
	delete[] rp, rg;
}

int main() {
		Initialize();
		ResetParameters();
		tbb::tick_count start = tbb::tick_count::now();
		tuple<double, double*>res = PSOSerial(TestFunction3, N);
		tbb::tick_count end = tbb::tick_count::now();
		cout << "Serial Time: " << (end - start).seconds() * 1000 << "ms" << endl;
		cout << "Min: " << get<0>(res) << endl;
		double* koord = get<1>(res);
		cout << "[ ";
		for (int i = 0; i < DIMENSIONS; i++)
		{
			cout << koord[i] << " ";
		}
		cout << "]" << endl;

		ResetParameters();

		start = tbb::tick_count::now();
		res = PSOParallel(TestFunction3, N);
		end = tbb::tick_count::now();
		cout << "Parallel Time: " << (end - start).seconds() * 1000 << "ms" << endl;
		cout << "Min: " << get<0>(res) << endl;
		koord = get<1>(res);
		cout << "[ ";
		for (int i = 0; i < DIMENSIONS; i++)
		{
			cout << koord[i] << " ";
		}
		cout << "]" << endl;
		Destruct();
}