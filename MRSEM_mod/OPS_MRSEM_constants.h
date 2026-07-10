double u0;
double simulation_time;
double dt;
double delta;
double r_max;
int ny;
int nz;
int niter;
int nregions;
double x_min;
double x_max;
double x_plane;
double y_min;
double y_max;

double z_min; 
double z_max;
double z_midpt;
double z_l; 
double eddy_z_min;
double eddy_z_max;
double Tinf;
double Twall;
double gama;
double Minf;

//double* x_radii;
double* eddy_vel;
double* radii;
int reg_num;
double* eddy_y_min;
double* eddy_y_max;
double* vols;

int* eddies;
int a;
int c;
int m;
int seed_gbl;

double* eddy_pos1 = nullptr; // array size 3 - 1 for each direction
double* eddy_r1 = nullptr; // array size 3 - 1 for each direction
int* eddy_eps1 = nullptr; // array size 3 - 1 for each direction

double* eddy_pos2 = nullptr;
double* eddy_r2 = nullptr;
int* eddy_eps2 = nullptr;

double* eddy_pos3 = nullptr;
double* eddy_r3 = nullptr;
int* eddy_eps3 = nullptr;

double* yinterp = nullptr;
double* r11interp = nullptr;
double* r21interp = nullptr;
double* r22interp = nullptr;
double* r33interp = nullptr;



