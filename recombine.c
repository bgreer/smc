#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"

void recombine ()
{
	int i, high;

	/* Find highest subdivision level */
	high = 0;
	for (i=0; i<N+1; i++)
		if (sd[i] > high) high = sd[i];

	for (i=high; i>=1; i/=2)
		recombine_level(i);
}

void recombine_level (int l)
{
	int i, j, count, done;
	/* sweep right */
	count= 0;
	done = 0;
	for (i=N-1; i>0; i--) /* search for first relevant point */
	{
		if (sd[i+1] == l/2 && sd[i] == l)
		{
			/* Look for points to the right that need recombination */
			count = 0;
			while (sd[i-count]==l && fabs(dx*dx*d2dx2(Tm,i-count)/(sd[i-count]*sd[i-count]))<1e-6)
				count++;
			i -= count;
			if (count>10) {done = 1;printf("count = %d, %d, %d\n", count, i, l);}
		}
	}
	if (done && 0) 
	{
		output_state("state", 0);
		exit(0);
	}
}
