double u0;
double dt;
double delta;
double r_max;
int ny;
int nz;
int niter;
double x_min;
double x_max;
double x_plane;
double y_min;
double y_max;
double eddy_y_min;
double eddy_y_max;
double z_min; 
double z_max;
double eddy_z_min;
double eddy_z_max;
double vol;
double rep_radius;
double dmin;
double gama;
double Minf;
double Twall;
double Tinf;
int eddies;
int a;
int c;
int m;
int seed_gbl;

// probability distribution variables
double powfact;
double cutoff;
double edgeval;
double vsf;
double rng_minval;
double rng_sf;
double k;
double b;

double* x_gbl = nullptr;
double* y_gbl = nullptr;
double* z_gbl = nullptr;
double* r_gbl = nullptr;
double* increment_gbl = nullptr;
int* eps_x_gbl = nullptr;
int* eps_y_gbl = nullptr;
int* eps_z_gbl = nullptr;

double* yinterp = nullptr;
double* r11interp = nullptr;
double* r21interp = nullptr;
double* r22interp = nullptr;
double* r33interp = nullptr;



