#include "solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

/** Retourne la différence (en secondes) entre deux timespec */
double get_delta(struct timespec begin, struct timespec end) {
	return end.tv_sec - begin.tv_sec + (end.tv_nsec - begin.tv_nsec) * 1e-9;
	}

int main(int argc, char * argv[]) {
	MPI_Init(&argc, &argv);
	int comm_rank, comm_size;
	struct timespec begin_steps, end_steps;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	if (argc < 3) {
		if (comm_rank == 0) {
			printf("USAGE: %s <nx> <ny>\n", argv[0]);
		}
		MPI_Finalize();
		return -1;
	}
	
	int nx = atoi(argv[1]);
	int ny = atoi(argv[2]);

	int NX = atoi(argv[3]); // Nombre de domaines selon X (Q3.2)
	int NY = comm_size / NX; // Nombre de domaines selon Y (Q1.3)
	int printMean = 0; //Fonction pour imprimer la moyenne à chaque étape

	heat_problem pb;
	create_problem(nx, ny, 1.0, NX, NY, &pb);

	double dt = 0.1*pb.dx*pb.dy;
	int niter = 1000;
	clock_gettime(CLOCK_REALTIME, &begin_steps);
	for (int i = 0; i < niter; i++) {
		step(&pb, dt);
		if (printMean && NX) { 
			print_mean(&pb);
		}
	}
	clock_gettime(CLOCK_REALTIME, &end_steps);
	
	if (comm_rank == 0) {
		printf("Temps de calcul: %lf s.\n", get_delta(begin_steps, end_steps));
	}
	
	print_result(&pb);

	free_problem(&pb);
	MPI_Finalize();
	return 0;
}

