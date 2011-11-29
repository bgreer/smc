#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"


void solve_pde ()
{
	int i;

	for (i=1; i<N; i++)
	{
		Tm[i] = Tm_old[i] + dt*d2dx2(Tm_old,i)
				- dt*ddx_2(W_old,Tf_old,i);
		Tf[i] = Tf_old[i] + dt*d2dx2(Tf_old,i)
				- dt*k*k*Tf_old[i] - dt*W_old[i]*ddx(Tm_old,i)
				- C*dt*W_old[i]*ddx(Tf_old,i)
				- C*dt*Tf_old[i]*ddx(W_old,i);
		Z[i] = Z_old[i] + p*dt*d2dx2(Z_old,i)
				- p*dt*k*k*Z_old[i] - p*dt*R*k*k*Tf_old[i]
				- C*2.0*dt*Z_old[i]*ddx(W_old,i)
				- C*dt*W_old[i]*ddx(Z_old,i);
	}
}

/* Find the first derivative of an array at a given point
	must take adaptive grid into account*/
double ddx (double* arr, int x0)
{
/*	return (arr[x0+1]-arr[x0-1])/(2.0*dx);*/
	if (x0==0 || x0==N) return 0.0;
	if (sd[x0] == sd[x0+1] && sd[x0] == sd[x0-1])
	{
		/* adjacent points are same resolution, use normal 1st deriv */
		return (arr[x0+1]-arr[x0-1])/(2.0*dx/sd[x0]);
	} else {
		/* on the edge of a regridding */
		if (sd[x0]/sd[x0-1]>2 || sd[x0-1]/sd[x0]>2 || sd[x0]/sd[x0+1]>2 || sd[x0+1]/sd[x0]>2)
		{
			printf("DDX ERROR: improper regridding at %d  (%d %d %d)\n", 
				x0, sd[x0-1], sd[x0], sd[x0+1]);
			exit(EXIT_FAILURE);
		}
		if (sd[x0] != sd[x0+1] && sd[x0] != sd[x0-1])
		{
			if (sd[x0+1]==sd[x0-1])
				return (arr[x0+1]-arr[x0-1])/(2.0*dx/sd[x0]);
			if (sd[x0-1] < sd[x0])
				return (arr[x0+1]-0.5*(arr[x0]+arr[x0-1]))/(2.0*dx/sd[x0]);
			else
				return (0.5*(arr[x0+1]+arr[x0])-arr[x0-1])/(2.0*dx/sd[x0]);
		}
		if (sd[x0]==sd[x0-1])
		{
			if (sd[x0+1] > sd[x0])
				return (arr[x0+1]-arr[x0-1])/(2.0*dx/sd[x0]);
			else
				return (0.5*(arr[x0+1]+arr[x0])-arr[x0-1])/(2.0*dx/sd[x0]);
		} else if (sd[x0]==sd[x0+1]) {
			if (sd[x0-1] > sd[x0])
				return (arr[x0+1]-arr[x0-1])/(2.0*dx/sd[x0]);
			else
				return (arr[x0+1]-0.5*(arr[x0-1]+arr[x0]))/(2.0*dx/sd[x0]);
		}
	}
}

double ddx_2 (double* arr1, double* arr2, int x0)
{
/*	return (arr[x0+1]-arr[x0-1])/(2.0*dx);*/
	if (x0==0 || x0==N) return 0.0;
	if (sd[x0] == sd[x0+1] && sd[x0] == sd[x0-1])
	{
		/* adjacent points are same resolution, use normal 1st deriv */
		return (arr1[x0+1]*arr2[x0+1]-arr1[x0-1]*arr2[x0-1])/(2.0*dx/sd[x0]);
	} else {
		/* on the edge of a regridding */
		if (sd[x0]/sd[x0-1]>2 || sd[x0-1]/sd[x0]>2 || sd[x0]/sd[x0+1]>2 || sd[x0+1]/sd[x0]>2)
		{
			printf("DDX_2 ERROR: improper regridding at %d  (%d %d %d)\n", 
				x0, sd[x0-1], sd[x0], sd[x0+1]);
			exit(EXIT_FAILURE);
		}
		if (sd[x0] != sd[x0+1] && sd[x0] != sd[x0-1])
		{
			if (sd[x0+1]==sd[x0-1])
				return (arr1[x0+1]*arr2[x0+1]-arr1[x0-1]*arr2[x0-1])/(2.0*dx/sd[x0]);
			if (sd[x0-1] < sd[x0])
				return (arr1[x0+1]*arr2[x0+1]-0.5*(arr1[x0]*arr2[x0]+arr1[x0-1]*arr2[x0-1]))/(2.0*dx/sd[x0]);
			else
				return (0.5*(arr1[x0+1]*arr2[x0+1]+arr1[x0]*arr2[x0])-arr1[x0-1]*arr2[x0-1])/(2.0*dx/sd[x0]);
		}
		if (sd[x0]==sd[x0-1])
		{
			if (sd[x0+1] > sd[x0])
				return (arr1[x0+1]*arr2[x0+1]-arr1[x0-1]*arr2[x0-1])/(2.0*dx/sd[x0]);
			else
				return (0.5*(arr1[x0+1]*arr2[x0+1]+arr1[x0]*arr2[x0])-arr1[x0-1]*arr2[x0-1])/(2.0*dx/sd[x0]);
		} else {
			if (sd[x0-1] > sd[x0])
				return (arr1[x0+1]*arr2[x0+1]-arr1[x0-1]*arr2[x0-1])/(2.0*dx/sd[x0]);
			else
				return (arr1[x0+1]*arr2[x0+1]-0.5*(arr1[x0-1]*arr2[x0-1]+arr1[x0]*arr2[x0]))/(2.0*dx/sd[x0]);
		}
	}
}

/* Find the second derivative of an array at a given point
	must take adaptive grid into account*/
double d2dx2 (double* arr, int x0)
{
/*	return (arr[x0+1]-2.0*arr[x0]+arr[x0-1])/(dx*dx);*/
	if (x0==0 || x0==N) return 0.0;
	if (sd[x0] == sd[x0+1] && sd[x0] == sd[x0-1])
	{
		/* adjacent points are same resolution, use normal 1st deriv */
		return (arr[x0+1]-2.0*arr[x0]+arr[x0-1])/(dx*dx/(sd[x0]*sd[x0]));
	} else {
		/* on the edge of a regridding */
		if (sd[x0]/sd[x0-1]>2 || sd[x0-1]/sd[x0]>2 || sd[x0]/sd[x0+1]>2 || sd[x0+1]/sd[x0]>2)
		{
			printf("D2DX2 ERROR: improper regridding at %d  (%d %d %d)\n", 
				x0, sd[x0-1], sd[x0], sd[x0+1]);
			exit(EXIT_FAILURE);
		}
		if (sd[x0] != sd[x0+1] && sd[x0] != sd[x0-1])
		{
			if (sd[x0-1]==sd[x0+1])
				return (arr[x0-1]-2.0*arr[x0]+arr[x0+1])/(dx*dx/(sd[x0]*sd[x0]));
			if (sd[x0-1] < sd[x0])
				return (0.5*(arr[x0-1]+arr[x0])-2.0*arr[x0]+arr[x0+1])/(dx*dx/(sd[x0]*sd[x0]));
			else
				return (arr[x0-1]-2.0*arr[x0]+0.5*(arr[x0+1]+arr[x0]))/(dx*dx/(sd[x0]*sd[x0]));
		}
		if (sd[x0]==sd[x0-1])
		{
			if (sd[x0+1] > sd[x0])
				return (arr[x0+1]-2.0*arr[x0]+arr[x0-1])/(dx*dx/(sd[x0]*sd[x0]));
			else
				return (0.5*(arr[x0+1]+arr[x0])-2.0*arr[x0]+arr[x0-1])/(dx*dx/(sd[x0]*sd[x0]));
		} else {
			if (sd[x0-1] > sd[x0])
				return (arr[x0+1]-2.0*arr[x0]+arr[x0-1])/(dx*dx/(sd[x0]*sd[x0]));
			else
				return (arr[x0+1]-2.0*arr[x0]+0.5*(arr[x0-1]+arr[x0]))/(dx*dx/(sd[x0]*sd[x0]));
		}
	}
}
