#include "solver.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void create_problem(int nx, int ny, double alpha, int NX, int NY, heat_problem * pb) {
	pb->nx = nx/NX + 2;
	pb->ny = ny/NY + 2;
	pb->dx = 1.0/(nx+1);
	pb->dy = 1.0/(ny+1);
	pb->alpha = alpha;
	pb->T = calloc(pb->nx*pb->ny, sizeof(double));
	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	pb->ycomm_rank = rank / NX;
	pb->ycomm_size = NY;
	pb->xcomm_rank = rank % NX;
	pb->xcomm_size = NX;
	
	if (pb->ycomm_rank == 0) {
		for (int j = 0; j < pb->nx; j++) {
			pb->T[j] = 0.0;
		}
	}

	if (pb->ycomm_rank == NY-1) {
		for (int j = 0; j < pb->nx; j++) {
			pb->T[(pb->ny-1)*pb->nx+j] = 0.0; 
		}
	}

	/*for (int i = 0; i < pb->ny; i++) {
		pb->T[i*pb->nx] = 1.0; 
		pb->T[i*pb->nx+(pb->nx-1)] = 1.0; 
	}
	*/
	
	if (pb->xcomm_rank == 0) {
		for (int i = 0; i < pb->ny; i++) {
			pb->T[i*pb->nx] = 1.0; 
		}
	}
	
	if (pb->xcomm_rank == NX - 1) {
		for (int i = 0; i < pb->ny; i++) {
			pb->T[i*pb->nx+(pb->nx-1)] = 1.0; 
		}
	}

}

void free_problem(heat_problem * pb) {
	free(pb->T);
}

void step(heat_problem * pb, double dt) {
	int nx = pb->nx;
	int ny = pb->ny;
	int size = nx * ny * sizeof(double);
	double * oldT = malloc(size);
	memcpy(oldT, pb->T, size);
	MPI_Request reqprev, reqnext;
	
	// Send pour les lignes
	if (pb->ycomm_rank != 0) {
		MPI_Isend(oldT + nx, nx, MPI_DOUBLE, (pb->ycomm_rank - 1)*pb->xcomm_size + pb->xcomm_rank, 0, MPI_COMM_WORLD, &reqprev);
		MPI_Recv(oldT, nx, MPI_DOUBLE, (pb->ycomm_rank - 1)*pb->xcomm_size + pb->xcomm_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Wait(&reqprev, MPI_STATUS_IGNORE);
	}
	if (pb->ycomm_rank != pb->ycomm_size-1) {
		MPI_Isend(oldT + nx*(ny - 2), nx, MPI_DOUBLE, (pb->ycomm_rank + 1)*pb->xcomm_size + pb->xcomm_rank, 0, MPI_COMM_WORLD, &reqnext);
		MPI_Recv(oldT + nx*(ny - 1), nx, MPI_DOUBLE, (pb->ycomm_rank + 1)*pb->xcomm_size + pb->xcomm_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Wait(&reqnext, MPI_STATUS_IGNORE);
	}
	
	// Send pour les colonnes
	if (pb->xcomm_rank != 0) {
		for (int i = 0; i < pb->ny; i++) {
			MPI_Isend(oldT + nx*i + 1, 1, MPI_DOUBLE, pb->ycomm_rank*pb->xcomm_size + pb->xcomm_rank - 1, 0, MPI_COMM_WORLD, &reqprev);
			MPI_Recv(oldT + nx*i, 1, MPI_DOUBLE, pb->ycomm_rank*pb->xcomm_size + pb->xcomm_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Wait(&reqprev, MPI_STATUS_IGNORE);
		}
	}
	if (pb->xcomm_rank != pb->xcomm_size-1) {	
		for (int i = 0; i < pb->ny; i++) {
			MPI_Isend(oldT + (nx - 2) + nx*i, 1, MPI_DOUBLE, pb->ycomm_rank*pb->xcomm_size + pb->xcomm_rank + 1, 0, MPI_COMM_WORLD, &reqnext);
			MPI_Recv(oldT + (nx - 1) + nx*i, 1, MPI_DOUBLE, pb->ycomm_rank*pb->xcomm_size + pb->xcomm_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Wait(&reqnext, MPI_STATUS_IGNORE);
		}
	}

	for (int i = 1; i < ny-1; i++) {
		for (int j = 1; j < nx-1; j++) {
			int index = i*nx + j;
			double lapl_x = (oldT[index+1] - 2*oldT[index] + oldT[index-1])/(pb->dx*pb->dx);
			double lapl_y = (oldT[index+nx] - 2*oldT[index] + oldT[index-nx])/(pb->dy*pb->dy);
			pb->T[index] += pb->alpha * dt * (lapl_x + lapl_y);
		}
	}
	free(oldT);
}

void print_result(heat_problem * pb) {

	if (pb->xcomm_size == 1) {
		int nx = (pb->nx - 2)*pb->xcomm_size + 2;
		int ny = (pb->ny - 2)*pb->ycomm_size + 2;
		int size = nx * ny * sizeof(double);
		double * newT = (double *)calloc(size, sizeof(double));

		MPI_Gather(pb->T + nx, nx*(pb->ny - 2), MPI_DOUBLE, newT + nx + nx*(pb->ny - 2)*pb->ycomm_rank, nx*(pb->ny - 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if (pb->ycomm_rank == 0) {
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++){
					printf("%3.2f ", newT[i*nx + j]);
				}
				printf("\n");
			}
		}
		
		free(newT);
	}
	
	else {
		/*
		int nx = (pb->nx - 2);
		int ny = (pb->ny - 2);
		int size = nx * pb->xcomm_size * ny * pb->ycomm_size * sizeof(double);
		double * newT = (double *)calloc(size, sizeof(double));
		
		for (int i=1; i <= ny; i++){
			MPI_Gather(pb->T + 1 + pb->nx*i, nx, MPI_DOUBLE, (newT + nx*pb->xcomm_size*ny*pb->ycomm_rank + pb->xcomm_rank*nx + nx*pb->xcomm_size*(i-1)), nx, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}
		*/

		MPI_Request* reqs = (MPI_Request*)malloc((pb->ny - 2) * sizeof(MPI_Request));
		
		for (int i=1; i < pb->ny - 1; i++){
			MPI_Isend(pb->T + 1 + pb->nx*i, pb->nx - 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &reqs[i-1]);
		}

		if (pb->ycomm_rank + pb->xcomm_rank == 0) {
			int nx = (pb->nx - 2)*pb->xcomm_size + 2;
			int ny = (pb->ny - 2)*pb->ycomm_size + 2;
			int size = nx * ny * sizeof(double);
			double * newT = (double *)calloc(size, sizeof(double));

			for (int r=0; r < pb->xcomm_size*pb->ycomm_size; r++){
				int ycomm_rank = r / pb->xcomm_size;
				int xcomm_rank = r % pb->xcomm_size;

				for (int i=1; i < pb->ny - 1; i++){
					MPI_Recv(newT + 1 + nx*(pb->ny - 2)*ycomm_rank + (pb->nx - 2)*xcomm_rank + nx*i, pb->nx - 2, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
				}
			}

			for (int i=1; i < ny - 1; i++) {
				newT[nx*i] = 1.0;
				newT[nx*(i + 1) - 1] = 1.0;
			}

			// Attendre que toutes les requêtes d'envoi soient terminées
			MPI_Waitall(pb->ny - 2, reqs, MPI_STATUSES_IGNORE);

			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < nx; j++){
					printf("%3.2f ", newT[i*nx + j]);
				}
				printf("\n");
			}

			free(newT);
		}

		free(reqs);
	}
}



void print_mean(heat_problem * pb) {
	if (pb->xcomm_size == 1) {
		int nx = pb->nx;
		int ny = pb->ny-2;
		double mean;
		double * MeanT = (double *)calloc(nx, sizeof(double));
		
		for (int j = 0; j < nx; j++) {
			mean = 0;
			for (int i = 0; i < ny; i++) {
				mean += pb->T[(i+1)*nx + j];
			}
			MeanT[j] = mean / ny;
		}
		
		double * NewMeanT = (double *)calloc(nx, sizeof(double));
		
		MPI_Reduce(MeanT, NewMeanT, nx, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		if (pb->ycomm_rank + pb->xcomm_rank == 0) {
			for (int j = 0; j < nx; j++){
				NewMeanT[j] = NewMeanT[j] / pb->ycomm_size;
				printf("%3.2f ", NewMeanT[j]);
			}
			printf("\n");
		}
		
		free(MeanT);
		free(NewMeanT);
	}
}

