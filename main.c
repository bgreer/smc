#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"

int main (int argc, char* argv[])
{
	int i, iter, subdivide_check, check_interval;
	int silent;
	double dt_base;

	/* Set default run params, overwrite with command line args */
	R = 10000.0;
	p = 1.0;
	k = 1.0;
	C = 0.0;

	N = 100;
	maxtime = 100.0;
	maxiter = 1000000;

	/* Parse command-line arguments */
	if (!parse_arguments(argc, argv))
	{
		usage();
		printf("Command-line argument parsing error. Quitting.\n");
		return EXIT_FAILURE;
	}

	dt_base = 0.01*p;
	maxsd = 1;
	dx = 1.0/N;
	dt = dt_base*dx*dx/(maxsd*maxsd); /* Attempt to keep things stable */

	time = 0.0;
	check_interval = 1000;
	debug2 = 0;

	silent = 1;

	/* Output run params */
	if (!silent)
	{
		printf("Using the following parameters:\n");
		printf("\tR = %6.1f\n", R);
		printf("\tp = %f\n", p);
		printf("\tk = %f\n", k);
		printf("\tC = %f\n", C);
		printf("\n\tGrid Resolution = %d\n", N);
		printf("\tTimestep = %e\n", dt);
	}

	/* Allocate and initialize arrays */
	init(initfile);
	boundary_cond();
	poisson_solve();

	/* Begin computation */
	iter = 0;
	subdivide_check = 0;
	do
	{
		if (iter==maxiter)
		{
			printf("Could not converge in %d iterations.\n", maxiter);
			break;
		}
		if (time >= maxtime)
		{
			printf("Max time reached: %f\n", time);
			break;
		}
		save_old();
/*		step_implicit();*/
		solve_pde();
		poisson_solve();
		boundary_cond();
		if (subdivide_check == check_interval)
		{
	/*		recombine();*/
			subdivide();
			subdivide_check = 0;
			dt = dt_base*dx*dx/(maxsd*maxsd);
		}
		iter++;
		if (debug2<0) break;
		if (debug1) debug2++;
		subdivide_check++;
		time += dt;
	} while (steady_state(0.0001)==0);
	if (!silent) printf("Solver took %d iterations to find solution (t=%f)\n", iter, time);

	if (!silent) printf("Final Nusselt number at top = %f\n", nusselt());

	if (silent) printf("%e\t%e\t%e\t%e\t%e\t%d\t%e\n", R, k, p, C, nusselt(), iter, time);

	output_state("state", 0);
	free(W);
	free(Z);
	free(Tf);
	free(Tm);
	return EXIT_SUCCESS;
}

void allocate_arrays ()
{
	sd = malloc((N+1)*sizeof(int));
	flag = malloc((N+1)*sizeof(int));

	W = malloc((N+1)*sizeof(double));
	Z = malloc((N+1)*sizeof(double));
	Tf = malloc((N+1)*sizeof(double));
	Tm = malloc((N+1)*sizeof(double));

	W_old = malloc((N+1)*sizeof(double));
	Z_old = malloc((N+1)*sizeof(double));
	Tf_old = malloc((N+1)*sizeof(double));
	Tm_old = malloc((N+1)*sizeof(double));
}

/*
	Set simple boundary conditions:
	- 0 velocity at walls
	- Bottom Tm set to 1.0
*/
void boundary_cond ()
{
	W[0] = W[N] = 0.0;
	//W[1] = W[N-1] = 0.0;
	Tm[0] = 1.0;
	Tm[N] = 0.0;
	Tf[0] = Tf[N] = 0.0;
	Z[0] = Z[N] = 0.0;
}

/*
	Copies arrays to _old arrays
	Do this before stepping in time
*/
void save_old ()
{
	int i;
	for (i=0; i<N+1; i++)
	{
		W_old[i] = W[i];
		Tf_old[i] = Tf[i];
		Tm_old[i] = Tm[i];
		Z_old[i] = Z[i];
	}
}

/*
	Initializes the 4 arrays to 0 (default)
	or to an input file in the format:

	N
	x	W[0]	Z[0]	Tf[0]	Tm[0]
	x	W[1]	Z[1]	Tf[1]	Tm[1]
	...
	x	W[N+1]	Z[N+1]	Tf[N+1]	Tm[N+1]
*/
void init (char* fname)
{
	int i, trash;
	double ftrash, ftrash2;
	float data1, data2, data3, data4, data5;
	FILE* fp;

	if (fname==NULL)
	{
		allocate_arrays();
		/* Set default initial conditions */
		for (i=0; i<N+1; i++)
		{
			Tf[i] = 0.0;
			Tm[i] = 1.0-(i*1.0)/(N+1.0);
			Z[i] = -0.01*sin(3.14*i/N);
			W[i] = 0.0;
			sd[i] = 1;
			flag[i] = 0;
		}
	} else {
		/* Load initial conditions from file */
		fp = fopen(fname, "r");
		i = 0;
		fscanf(fp, "%d\n", &trash);
		N = trash;
		allocate_arrays();
		for (i=0; i<N+1; i++)
		{
			fscanf(fp, "%e\t%e\t%e\t%e\t%e\t%e\t%d\n", 
				&ftrash, &data1, &data2, &data3, &data4, &ftrash2, &trash);
			W[i] = data1;
			Z[i] = data2;
			Tf[i] = data3;
			Tm[i] = data4;
			sd[i] = trash;
			flag[i] = 0;
		}
		fclose(fp);
		printf("Init successfully read from %s\n", fname);
	}
}

void output_state (char* fname, int version)
{
	char* buf;
	FILE* fp;
	int i;
	double pos;

	if (fname==NULL)
	{
		printf("ERROR: Must supply filename for state output\n");
		return;
	}

	buf = malloc(50*sizeof(char));
	sprintf(buf, "%s%.5d", fname, version);

	fp = fopen(buf, "w");
	fprintf(fp, "%d\n", N);
	pos = 0.0;
	for (i=0; i<N+1; i++)
	{
		fprintf(fp, "%e\t%e\t%e\t%e\t%e\t%e\t%d\n", pos, W[i], Z[i], Tf[i], Tm[i]
			, d2dx2(Tm,i), sd[i]);
		if (i<N)
		{
			if (sd[i+1] > sd[i])
				pos += dx/(1.0*sd[i]);
			else
				pos += dx/(1.0*sd[i+1]);
		}
	}
	fclose(fp);
	free(buf);
}

/*
	Used to parse command-line arguments
	User can specify the following:
	-R #    Rayleigh number
	-p #    Prandtl number
	-k #    Horizontal wavenumber
	-C #    Self-interaction parameter
*/
int parse_arguments (int argc, char* argv[])
{
	int i;
	if (argc==1) return 1;
	if (argc%2==0) return 0;

	initfile = NULL;

	for (i=1; i<argc; i++)
	{
		if (strcmp(argv[i], "-R")==0)
		{
			R = atof(argv[i+1]);
			i++;
		} else if (strcmp(argv[i], "-p")==0) {
			p = atof(argv[i+1]);
			i++;
		} else if (strcmp(argv[i], "-k")==0) {
			k = atof(argv[i+1]);
			i++;
		} else if (strcmp(argv[i], "-C")==0) {
			C = atof(argv[i+1]);
			i++;
		} else if (strcmp(argv[i], "-init")==0) {
			initfile = argv[i+1];
			i++;
		} else if (strcmp(argv[i], "-maxiter")==0) {
			maxiter = atoi(argv[i+1]);
			i++;
		} else if (strcmp(argv[i], "-maxtime")==0) {
			maxtime = atof(argv[i+1]);
			i++;
		} else 
			return 0;
	}
}

void usage ()
{
	printf("Incorrect Usage.\n");
}
