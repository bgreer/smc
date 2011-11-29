int debug1, debug2;

/* Numerical Properties */
int N; /* Grid Size */
double dx, dt; /* Base spatial and temporal steps */
int *sd, *sd2; /* List of resolution level per point */
int *flag, *flag2, *newpt; /* Flag for re-gridding */
int maxsd;

/* Constants */
double R; /* Rayleigh number */
double p; /* Prandtl number */
double k; /* Horizontal wavenumber */
double C; /* Self-interaction parameter */

/* Physical Quantity Arrays */
double *W, *W_old; /* Vertical velocity */
double *Tf, *Tf_old; /* Fluctuating Temp */
double *Tm, *Tm_old; /* Mean Temp */
double *Z, *Z_old; /* Laplacian of W */

double time, maxtime;
int maxiter;
char* initfile;

/* Prototypes */
void boundary_cond ();
void save_old ();
void init (char* fname);
void output_state (char* name, int version);
int parse_arguments (int argc, char* argv[]);
void usage ();

/* poisson.c */
void poisson_solve ();
void poisson_solve_old ();

/* solver.c */
void solve_pde ();
double ddx (double* arr, int x0);
double ddx_2 (double* arr1, double* arr2, int x0);
double d2dx2 (double* arr, int x0);

/* implicit.c */
void step_implicit ();
void matrix_add (double *A, double *B, double *X);

/* stats.c */
double nusselt ();
int steady_state (double tol);
void subdivide ();
void spread_flags (int i);
void subdivide_level (int l);
void expand (int l);
void stuff ();

/* recombine.c */
void recombine ();
void recombine_level (int l);
