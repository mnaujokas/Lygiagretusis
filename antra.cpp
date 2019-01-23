#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>

using namespace std;

int numDP = 1000;      // Vietoviu skaicius (demand points, max 10000)
int numPF = 10;        // Esanciu objektu skaicius (preexisting facilities)
int numCL = 24;        // Kandidatu naujiems objektams skaicius (candidate locations)
int numX;         // Nauju objektu skaicius
int world_rank;
int world_size;
int num_elements_per_proc;
double **demandPoints; // Geografiniai duomenys


//=============================================================================

double getTime();
void loadDemandPoints();
double HaversineDistance(double* a, double* b);
double evaluateSolution(int *X);
int increaseX(int *X, int index, int maxindex);
void calculate();

//=============================================================================

int main(){
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	num_elements_per_proc = numDP/world_size;
 	double ts = getTime();
	loadDemandPoints();

 	numX = 3;
 	calculate();
 	numX = 4;
 	calculate();
 	numX = 5;
 	calculate();

 	double tf = getTime();
 	if(world_rank==0)cout<<"Visu skaiciavimu laikas vienu metu: " << tf - ts << endl;
	MPI_Finalize();
}


void calculate() {

	double ts = getTime();          // Algoritmo vykdymo pradzios laikas

	          // Nuskaitomi duomenys
	// Sudarom pradini sprendini: [0, 1, 2, 3, ...]
	int *X = new int[numX];
	int *bestX = new int[numX];
	for (int i=0; i<numX; i++) {
		X[i] = i;
		bestX[i] = i;


	}
	//Sarom iteracijusarasa
	int *data = new int[numDP];

	if(world_rank==0){
		for (int i=0; i<numDP; i++) {
			data[i] = i;
		}
	}
	//Sudarom resultatu sarasa
	double *results = new double[world_size];
	//sudarom iteraciju poaibi
	int *sub_data = new int[num_elements_per_proc];

	double u = evaluateSolution(X);
	double bestU = u;

	//----- Pagrindinis ciklas ------------------------------------------------
	//padalinam iteracijas
	MPI_Scatter(data, num_elements_per_proc, MPI_INT, sub_data, num_elements_per_proc, MPI_INT, 0, MPI_COMM_WORLD);
	while (true) {
		if (increaseX(X, numX-1, numCL)) {
			double U = 0;
			int bestPF;
			int bestXE;
			double d;
			for (int i=0; i<num_elements_per_proc; i++) {
				bestPF = 1e5;
				for (int j=0; j<numPF; j++) {
					d = HaversineDistance(demandPoints[sub_data[i]], demandPoints[j]);
					if (d < bestPF) bestPF = d;
				}
				bestXE = 1e5;
				for (int j=0; j<numX; j++) {
					d = HaversineDistance(demandPoints[sub_data[i]], demandPoints[X[j]]);
					if (d < bestXE) bestXE = d;
				}
				if (bestXE < bestPF) U += demandPoints[sub_data[i]][2];
				else if (bestXE == bestPF) U += 0.3*demandPoints[sub_data[i]][2];
			}

			MPI_Gather(&U, 1, MPI_DOUBLE, results, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			//u=U;
			if(world_rank==0){
				for (int i=0; i<world_size; i++){
					if (u<=results[i]) u=results[i];
				}

				if (u > bestU) {
					bestU = u;
					for (int i=0; i<numX; i++) bestX[i] = X[i];
				}
			}

		}
		else break;
	}
	//----- Rezultatu spausdinimas --------------------------------------------

	;     // Skaiciavimu pabaigos laikas
	if(world_rank==0){
		double tf = getTime();
		cout << "Geriausias sprendinys: ";
		for (int i=0; i<numX; i++) cout << bestX[i] << " ";
		cout << "(" << bestU << ")" << endl;
		cout << "Skaiciavimo trukme: " << tf-ts << endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);

}

//=============================================================================

void loadDemandPoints() {
	FILE *f;
	f = fopen("demandPoints.dat", "r");
	demandPoints = new double*[numDP];
	for (int i=0; i<numDP; i++) {
		demandPoints[i] = new double[3];
		fscanf(f, "%lf%lf%lf", &demandPoints[i][0], &demandPoints[i][1], &demandPoints[i][2]);
	}
	fclose(f);
}

//=============================================================================

double HaversineDistance(double* a, double* b) {
   double dlon = fabs(a[0] - b[0]);
   double dlat = fabs(a[1] - b[1]);
   double aa = pow((sin((double)dlon/(double)2*0.01745)),2) + cos(a[0]*0.01745) * cos(b[0]*0.01745) * pow((sin((double)dlat/(double)2*0.01745)),2);
   double c = 2 * atan2(sqrt(aa), sqrt(1-aa));
   double d = 6371 * c;
   return d;
}

//=============================================================================

double getTime() {
   struct timeval laikas;
   gettimeofday(&laikas, NULL);
   double rez = (double)laikas.tv_sec+(double)laikas.tv_usec/1000000;
   return rez;
}

//=============================================================================

double evaluateSolution(int *X) {
	double U = 0;
	int bestPF;
	int bestX;
	double d;
	for (int i=0; i<numDP; i++) {
		bestPF = 1e5;
		for (int j=0; j<numPF; j++) {
			d = HaversineDistance(demandPoints[i], demandPoints[j]);
			if (d < bestPF) bestPF = d;
		}
		bestX = 1e5;
		for (int j=0; j<numX; j++) {
			d = HaversineDistance(demandPoints[i], demandPoints[X[j]]);
			if (d < bestX) bestX = d;
		}
		if (bestX < bestPF) U += demandPoints[i][2];
		else if (bestX == bestPF) U += 0.3*demandPoints[i][2];
	}
	return U;
}

//=============================================================================

int increaseX(int *X, int index, int maxindex) {
	if (X[index]+1 < maxindex-(numX-index-1)) {
		X[index]++;
	}
	else {
		if ((index == 0) && (X[index]+1 == maxindex-(numX-index-1))) {
			return 0;
		}
		else {
			if (increaseX(X, index-1, maxindex)) X[index] = X[index-1]+1;
			else return 0;
		}
	}
	return 1;
}
