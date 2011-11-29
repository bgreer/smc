#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"

/*
	Computes the Nusselt number at the upper boundary
*/
double nusselt ()
{
	return -ddx(Tm,N-1) + W[N-1]*Tf[N-1];
}

/*
	Determines if the solutions have reached steady state

*/
int steady_state (double tol)
{
	int i;
	for (i=0; i<N+1; i++)
	{
		//printf("%d\tabs=%e\ttol=%e\n", i, fabs(Tm_old[i]-Tm[i]), tol*dt);
		if (fabs(W_old[i]-W[i]) > tol*dt) return 0;
		if (fabs(Tm_old[i]-Tm[i]) > tol*dt) return 0;
		if (fabs(Tf_old[i]-Tf[i]) > tol*dt) return 0;
	}
	return 1;
}

/*
	Determine if grid needs to be subdivided, if so, do it
*/
void subdivide ()
{
	int i, end = 0;
	int newcount, high;

	if (N >= 400) return;

	/* Scan through to find regions that need refinement */
	newcount = 0;
	for (i=1; i<N; i++)
	{
		if (fabs(dx*dx*d2dx2(Tm,i)/(sd[i]*sd[i])) > 0.01 && sd[i] <= 2)
		{
			/* subdivide between i-1, i, i+1 */
			/* printf("subdivide at %d, %f\n", i, fabs(dx*dx*d2dx2(Tm,i))/(sd[i]*sd[i]));*/
			flag[i] = 1;
			newcount++;
		}
	}
	if (newcount == 0) return;

	/* Spread regridding regions to prevent hard jumps */
	for (i=1; i<N; i++)
	{
		spread_flags(i);
	}

	/* find highest regridding level */
	high = 0;
	for (i=1; i<N; i++)
		if (flag[i] && sd[i] > high) high = sd[i];
	maxsd = high*2;

	/* Free memory for creating new things */
	free(Tm_old);
	free(W_old);
	free(Tf_old);
	free(Z_old);

	/* Subdivide starting with lowest level */
	for (i=1; i<=high; i*=2)
		subdivide_level(i);

	/* re-allocate the _old arrrays for further computation */
	Tm_old = malloc((N+1)*sizeof(double));
	Tf_old = malloc((N+1)*sizeof(double));
	W_old = malloc((N+1)*sizeof(double));
	Z_old = malloc((N+1)*sizeof(double));

	if (end)
	{
		stuff();
		output_state("state", 0);
		exit(EXIT_FAILURE);
	}
	if (N>100) debug1 = 1;

	/* clean up flags */
	for (i=0; i<N+1; i++)
		flag[i] = 0;
}

/* Spread flags to maintain proper gridding, called recursively */
void spread_flags (int i)
{
	int cont;
	if (!flag[i] || i==0 || i==N) return;
	/*printf("spreading at %d --(%d %d)  -(%d %d)  (%d %d)  +(%d %d) ++(%d %d)\n", i, flag[i-2], sd[i-2], flag[i-1], sd[i-1], flag[i], sd[i], flag[i+1], sd[i+1], flag[i+2], sd[i+2]);*/
	/* Left */
	cont = 1;
	if (sd[i] > sd[i-1] && !flag[i-1])
	{
		flag[i-1] = 1;
		spread_flags(i-1);
		cont = 0;
	} else if (cont && sd[i] == sd[i-1]) {
		if (i > 1)
		{
			if (sd[i] > sd[i-2] && !flag[i-2])
			{
				flag[i-2] = 1;
				spread_flags(i-2);
			} else if (sd[i] < sd[i-2] && !flag[i-1]) {
				flag[i-1] = 1;
				spread_flags(i-1);
			}
		}
	}

	/* Right */
	cont = 1;
	if (sd[i] > sd[i+1] && !flag[i+1])
	{
		flag[i+1] = 1;
		spread_flags(i+1);
		cont = 0;
	} else if (cont && sd[i] == sd[i+1]) {
		if (i < N-1)
		{
			if (sd[i] > sd[i+2] && !flag[i+2])
			{
				flag[i+2] = 1;
				spread_flags(i+2);
			} else if (sd[i] < sd[i+2] && !flag[i+1]) {
				flag[i+1] = 1;
				spread_flags(i+1);
			}
		}
	}
}

void subdivide_level (int l)
{
	int i, j, newN;
/*	printf("Subdivide at level %d\n", l);*/
	newN = N;

	if (flag[0] && sd[0]==l)
		newN++;
	if (flag[N] && sd[N]==l)
		newN++;
	for (i=1; i<N; i++)
	{
		/* Search for start of block */
		if (flag[i] && sd[i]==l)
		{
			newN += 2;
			/* Keep scanning for contiguous block */
			i++;
			while (flag[i] && sd[i]==l)
			{
				i++;
				newN++;
			}
		}
	}
/*	printf("Old N: %d  New N: %d\n", N, newN);*/

	if (newN==N)
	{
/*		printf("No need to subdivide\n");*/
		return;
	}

	/* Allocate Memory */
	flag2 = malloc((newN+1)*sizeof(int));
	sd2 = malloc((newN+1)*sizeof(int));
	newpt = malloc((newN+1)*sizeof(int));
	Tm_old = malloc((newN+1)*sizeof(double));
	Z_old = malloc((newN+1)*sizeof(double));
	Tf_old = malloc((newN+1)*sizeof(double));
	W_old = malloc((newN+1)*sizeof(double));

	/* Begin mapping old grid onto new grid (interpolation) */
	j = 0;

	for (i=0; i<N+1; i++)
	{
		if (flag[i] && sd[i]==l)
		{
			if (i>0)
			{
				Tm_old[j] = 0.5*(Tm[i]+Tm[i-1]);
				Tf_old[j] = 0.5*(Tf[i]+Tf[i-1]);
				W_old[j] = 0.5*(W[i]+W[i-1]);
				Z_old[j] = 0.5*(Z[i]+Z[i-1]);
				flag2[j] = 0;
				sd2[j] = l*2;
				newpt[j] = 1;
				j++;
			}
			Tm_old[j] = Tm[i];
			Tf_old[j] = Tf[i];
			W_old[j] = W[i];
			Z_old[j] = Z[i];
			flag2[j] = 0;
			sd2[j] = l*2;
			newpt[j] = 0;
			if (i<N)
			{
				if (!flag[i+1] || sd[i+1] != l)
				{
					j++;
					Tm_old[j] = 0.5*(Tm[i]+Tm[i+1]);
					Tf_old[j] = 0.5*(Tf[i]+Tf[i+1]);
					W_old[j] = 0.5*(W[i]+W[i+1]);
					Z_old[j] = 0.5*(Z[i]+Z[i+1]);
					flag2[j] = 0;
					sd2[j] = l*2;
					newpt[j] = 1;
				}
			}
		} else {
			Tm_old[j] = Tm[i];
			Tf_old[j] = Tf[i];
			W_old[j] = W[i];
			Z_old[j] = Z[i];
			flag2[j] = flag[i];
			sd2[j] = sd[i];
			newpt[j] = 0;
		}
		j++;
	}
/*	printf("remapping done, time to swap %d  %d\n", j, newN);*/
	N = newN;
/*	free(sd);*/
/*	free(flag);*/
	sd = sd2;
	flag = flag2;
/*	printf("swap complete\n");*/
	/* copy _old arrays to normal ones */
	Tm = Tm_old;
	Tf = Tf_old;
	W = W_old;
	Z = Z_old;
	expand(l);
/*	printf("Done subdividing\n");*/
}

/* Expand subdivision markers after a regridding */
void expand (int l)
{
	int i, *new;
	new = malloc((N+1)*sizeof(int));

	if (newpt[1]) new[0] = sd[1];
		else new[0] = sd[0];
	if (newpt[N-1]) new[N] = sd[N-1];
		else new[N] = sd[N];

	for (i=1; i<N; i++)
	{
		if ((sd[i] == sd[i-1]/2. && newpt[i-1]) || (sd[i] == sd[i+1]/2. && newpt[i+1]))
			new[i] = sd[i]*2;
		else
			new[i] = sd[i];
	}
	for (i=0; i<=N; i++)
		sd[i] = new[i];
	free(new);
}

void stuff()
{
	int i;
	for (i=0; i<=N; i++)
	{
/*		printf("%d\t%d\t%d\n", i, flag[i], sd[i]);*/
	}
}
