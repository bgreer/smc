#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_cblas.h>
#include "header.h"


void step_implicit ()
{
	int i;
	double fn, fp, dxs;
	double **A, **B, **V, **D, **x, *compresult, *compresult2;
	double *Tm_pr, *Tf_pr, *W_pr, *Z_pr;

	Tm_pr = malloc((N+1)*sizeof(double));
	Tf_pr = malloc((N+1)*sizeof(double));
	W_pr = malloc((N+1)*sizeof(double));
	Z_pr = malloc((N+1)*sizeof(double));
	memcpy(Tm_pr, Tm, (N+1)*sizeof(double));
	memcpy(Tf_pr, Tf, (N+1)*sizeof(double));
	memcpy(W_pr, W, (N+1)*sizeof(double));
	memcpy(Z_pr, Z, (N+1)*sizeof(double));

	/* Allocate required memory */
	A = malloc((N-2)*sizeof(double*));
	B = malloc((N-1)*sizeof(double*));
	V = malloc((N-2)*sizeof(double*));
	D = malloc((N-1)*sizeof(double*));
	x = malloc((N-1)*sizeof(double*));
	for (i=0; i<N-2; i++)
	{
		A[i] = calloc(16, sizeof(double));
		B[i] = calloc(16, sizeof(double));
		V[i] = calloc(16, sizeof(double));
		D[i] = calloc(4, sizeof(double));
		x[i] = calloc(4, sizeof(double));
	}
	B[N-2] = calloc(16, sizeof(double));
	D[N-2] = calloc(4, sizeof(double));
	x[N-2] = calloc(4, sizeof(double));

	compresult = calloc(16, sizeof(double));
	compresult2 = calloc(16, sizeof(double));

	/* Fill D array */
	for (i=0; i<N-1; i++)
	{
		D[i][0] = 0.5*d2dx2(Tm_old,i+1)-0.5*ddx_2(W_old,Tf_old,i+1)+Tm_old[i+1]/dt;
		D[i][1] = 0.0;
		D[i][2] = 0.5*(d2dx2(Tf_old,i+1) - k*k*Tf_old[i+1] - W_old[i+1]*ddx(Tm_old,i+1))
					- C*(W_old[i+1]*ddx(Tf_old,i+1) + 0.5*Tf_old[i+1]*ddx(W_old,i+1)) + Tf_old[i+1]/dt;
		D[i][3] = 0.5*p*(d2dx2(Z_old,i+1)-k*k*Z_old[i+1]-R*k*k*Tf_old[i+1])
					- C*(Z_old[i+1]*ddx(W_old,i+1)+0.5*W_old[i+1]*ddx(Z_old,i+1)) + Z_old[i+1]/dt;
	}

	/* Create Jacobian */
	for (i=0; i<N-2; i++)
	{
		if (sd[i+1] < sd[i+2])
			fn = 0.5;
		else
			fn = 1.0;
		if (sd[i+3] < sd[i+2])
			fp = 0.5;
		else
			fp = 1.0;
		dxs = dx/sd[i+2];

		A[i][0] = fn/(2.*dxs*dxs);
		A[i][1] = Tf_pr[i+1]*fn/(4.*dxs);
		A[i][2] = W_pr[i+1]*fn/(4.*dxs);
		A[i][3] = 0.0;

		A[i][4] = 0.0;
		A[i][5] = fn/(dxs*dxs);
		A[i][6] = 0.0;
		A[i][7] = 0.0;

		A[i][8] = W_pr[i+2]*fn/(4.*dxs);
		A[i][9] = Tf_pr[i+2]*fp*C/(4.*dxs);
		A[i][10] = fn/(2.*dxs*dxs)+W_pr[i+2]*fn*C/(2.*dxs);
		A[i][11] = 0.0;

		A[i][12] = 0.0;
		A[i][13] = Z_pr[i+2]*C*fn/(2.*dxs);
		A[i][14] = 0.0;
		A[i][15] = fn*p/(2.*dxs*dxs)+W_pr[i+2]*fn*C/(4.*dxs);
	}
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
		dxs = dx/sd[i+1];

		B[i][0] = -((fp+fn)/(2.*dxs*dxs)+1./dt);
		B[i][1] = -Tf_pr[i+1]*(fp-fn)/(4.*dxs);
		B[i][2] = -W_pr[i+1]*(fp-fn)/(4.*dxs);
		B[i][3] = 0.0;

		B[i][4] = 0.0;
		B[i][5] = -((fp+fn)/(dxs*dxs)+k*k);
		B[i][6] = 0.0;
		B[i][7] = -1.0;

		B[i][8] = -W_pr[i+1]*(fp-fn)/(4.*dxs);
		B[i][9] = Tf_pr[i]*fn*C/(2.*dxs)-Tf_pr[i+1]*C*(fp-fn)/(4.*dxs)-Tf_pr[i+2]*fp*C/(2.*dxs)+Tm_pr[i]*fn/(4.*dxs)
					- Tm_pr[i+1]*(fp-fn)/(4.*dxs) - Tm_pr[i+2]*fp/(4.*dxs);
		B[i][10] = -(fp+fn)/(2.*dxs*dxs)-k*k-1./dt-W_pr[i+1]*C*(fp-fn)/(4.*dxs)+W_pr[i]*fn*C/(4.*dxs)
					- W_pr[i+2]*fp*C/(4.*dxs);
		B[i][11] = 0.0;

		B[i][12] = 0.0;
		B[i][13] = Z_pr[i]*fn*C/(4.*dxs)-Z_pr[i+1]*3*C*(fp-fn)/(4.*dxs)-Z_pr[i+2]*C*fp/(4.*dxs);
		B[i][14] = -0.5*p*R*k*k;
		B[i][15] = -p*(fp+fn)/(2.*dxs*dxs)-1./dt-p*k*k/2.-W_pr[i+1]*3.*C*(fp-fn)/(4.*dxs)+W_pr[i]*C*fn/(2.*dxs)
					-W_pr[i+2]*C*fp/(2.*dxs);
	}
	for (i=0; i<N-2; i++)
	{
		if (sd[i] < sd[i+1])
			fn = 0.5;
		else
			fn = 1.0;
		if (sd[i+2] < sd[i+1])
			fp = 0.5;
		else
			fp = 1.0;
		dxs = dx/sd[i+1];

		V[i][0] = fp/(2.*dxs*dxs);
		V[i][1] = -Tf_pr[i+2]*fp/(4.*dxs);
		V[i][2] = -W_pr[i+2]*fp/(4.*dxs);
		V[i][3] = 0.0;

		V[i][4] = 0.0;
		V[i][5] = fp/(dxs*dxs);
		V[i][6] = 0.0;
		V[i][7] = 0.0;

		V[i][8] = -W_pr[i+1]*fp/(4.*dxs);
		V[i][9] = -Tf_pr[i+1]*fn*C/(4.*dxs);
		V[i][10] = fp/(2.*dxs*dxs)-W_pr[i+1]*C*fp/(2.*dxs);
		V[i][11] = 0.0;
		
		V[i][12] = 0.0;
		V[i][13] = -Z_pr[i+1]*C*fp/(2.*dxs);
		V[i][14] = 0.0;
		V[i][15] = p*fp/(2.*dxs*dxs)-W_pr[i+1]*C*fp/(4.*dxs);
	}

	/* Begin first pass */
	for (i=1; i<N-1; i++)
	{
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4,4,4,1.0, B[i], 4, B[i-1], 4, 0.0, compresult, 4);
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4,4,4,-1.0, A[i], 4, V[i-1], 4, 0.0, compresult2, 4);
		matrix_add(compresult, compresult2, B[i]);
		
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4,4,4,1.0, V[i], 4, B[i-1], 4, 0.0, compresult, 4);
		memcpy(V[i], compresult, 16*sizeof(double));

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4,4,4,1.0, D[i], 4, B[i-1], 4, 0.0, compresult, 4);
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4,4,4,-1.0, A[i], 4, D[i-1], 4, 0.0, compresult2, 4);
		matrix_add(compresult, compresult2, D[i]);
	}

	/* Back-sub */
	

	/* Free memory */
	free(A);free(B);free(V);free(D);free(x);
	free(Tm_pr);free(Tf_pr);free(W_pr);free(Z_pr);
}

/* COmputes X = A+B */
void matrix_add (double *A, double *B, double *X)
{
	int i;
	for (i=0; i<16; i++)
		X[i] = A[i]+B[i];
}
