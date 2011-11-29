#include <stdlib.h>
#include "header.h"

/*
	Given an array Z, finds the array W based on
	Z = d2/dx2[W] - k^2*W

	Uses the globally defined Z, W, k

	This method uses central difference for d2/dx2
*/
void poisson_solve ()
{
	int i;
	double *a, *b, *c, *d;
	double *cp, *dp, fp, fn;

	a = malloc((N-2)*sizeof(double));
	b = malloc((N-1)*sizeof(double));
	c = malloc((N-2)*sizeof(double));
	d = malloc((N-1)*sizeof(double));

	cp = malloc((N-2)*sizeof(double));
	dp = malloc((N-1)*sizeof(double));

	fp = fn = 1.0;

	/* Initialize a,b,c,d arrays */

	/* a */
	for (i=0; i<N-2; i++)
	{
		if (sd[i+1] < sd[i+2])
			fn = 0.5;
		else
			fn = 1.0;
		a[i] = fn;
	}
	/* b */
	for (i=0; i<N-1; i++)
	{
		if (sd[i] < sd[i+1])
			fn = 0.5;
		else
			fn = 1.0;
		if (sd[i+2] < sd[i+1])
			fp = 0.5;
		else
			fp = 1.0;
		b[i] = -(fp+fn+k*k*dx*dx/(sd[i+1]*sd[i+1]));
	}
	/* c */
	for (i=0; i<N-2; i++)
	{
		if (sd[i+2] < sd[i+1])
			fp = 0.5;
		else
			fp = 1.0;
		c[i] = fp;
	}
	/* d */
	for (i=0; i<N-1; i++)
	{
		d[i] = Z[i+1]*dx*dx/(sd[i+1]*sd[i+1]);
	}

	/* Create modified coeffs */
	cp[0] = c[0]/b[0];
	for (i=1; i<N-2; i++)
	{
		cp[i] = c[i]/(b[i]-cp[i-1]*a[i-1]);
	}

	dp[0] = d[0]/b[0];
	for (i=1; i<N-1; i++)
	{
		dp[i] = (d[i]-dp[i-1]*a[i-1])/(b[i]-cp[i-1]*a[i-1]);
	}

	/* Begin back-sub */
	W[N-1] = dp[N-2];
	for (i=N-3; i>=0; i--)
	{
		W[i+1] = dp[i] - cp[i]*W[i+2];
	}

	/* Free memory like a good programmer */
	free(a);
	free(b);
	free(c);
	free(d);
	free(cp);
	free(dp);
}

/* For testing purposes only */
void poisson_solve_old ()
{
	int i;
	float *c, *d;

	c = malloc((N-2)*sizeof(float));
	d = malloc((N-1)*sizeof(float));

	/* Create modified coeffs */
	c[0] = -1.0/(2.0+k*k*dx*dx);
	for (i=1; i<N-2; i++)
	{
		c[i] = -1.0/((2.0+k*k*dx*dx)+c[i-1]);
	}

	d[0] = -Z[1]*dx*dx/(2.0+k*k*dx*dx);
	for (i=1; i<N-1; i++)
	{
		d[i] = -(Z[i+1]*dx*dx-d[i-1])/((2.0+k*k*dx*dx)+c[i-1]);
	}

	/* Begin back-sub */
	W[N-1] = d[N-2];
	for (i=N-3; i>=0; i--)
	{
		W[i+1] = d[i] - c[i]*W[i+2];
	}

	/* Free memory like a good programmer */
	free(c);
	free(d);

}
