smc : main.o solver.o poisson.o stats.o recombine.o implicit.o header.h
	gcc -O2 main.o solver.o poisson.o stats.o recombine.o implicit.o -o smc -lm -lgslcblas

main.o : main.c header.h
	gcc -c main.c header.h

solver.o : solver.c header.h
	gcc -c -O2 solver.c header.h

poisson.o : poisson.c header.h
	gcc -c -O2 poisson.c header.h

stats.o : stats.c header.h
	gcc -c -O2 stats.c header.h

recombine.o : recombine.c header.h
	gcc -c -O2 recombine.c header.h

implicit.o : implicit.c header.h
	gcc -c implicit.c header.h
