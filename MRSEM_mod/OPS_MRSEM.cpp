#include <cstdlib>
#include <cmath>
#include <random>
#include <map>
#include <stdio.h>
#include <fstream>
#include <string>


#define OPS_2D
#include "ops_seq.h"
#include "ops_lib_core.h"
#include "TBL_data.h"
#include "OPS_MRSEM_constants.h"
#include "OPS_MRSEM_eddy_functions.h"
#include "OPS_MRSEM_kernels.h"


int main(int argc, char** argv){

    ops_init(argc, argv, 1);

    // ------------------------------- Declare constant variable values ------------------------------------------
    u0 = 823.6;
    simulation_time = 0.0;
    dt = 0.00000002 * 5.0;
    delta = 0.007;
    r_max = 0.41 * delta;
    ny = 100;
    nz = 150;
    niter = 2000;
    nregions = 3;
    x_min = -r_max;
    x_max = r_max;
    x_plane = 0.0;
    y_min = 0.0;
    y_max = 0.007; //1.5 * delta + r_max;

    z_min = 0.0;
    z_max = 0.05;
    z_midpt = z_min + 0.5 * (z_max - z_min);
    z_l = z_max - z_min;
    eddy_z_min = z_min - r_max;
    eddy_z_max = z_max + r_max;
    Tinf = 202.2;
    Twall = 1.7 * Tinf;
    gama = 1.4;
    Minf = 2.0;

    eddies = (int*)malloc(nregions * sizeof(int));

    eddy_vel = (double*)malloc(sizeof(double) * nregions);
    eddy_vel[0] = 0.62 * u0;
    eddy_vel[1] = 0.62 * u0;
    eddy_vel[2] = 0.8 * u0;

    radii = (double*)malloc(nregions * 3 * sizeof(double));
    reg_num = 0;

    


    // each 3 set is t, y, z of the ith region
    // initially assign x lengths for eddy calculation, then switch to time after number of eddies has been calculated
    double lplus = 3.4178e-5;
    
    radii[0] = 100.0 * lplus;
    radii[1] = 20.0 * lplus;
    radii[2] = 60.0 * lplus;
    
    radii[3] = 60.0 * lplus;
    radii[4] = 30.0 * lplus;
    radii[5] = 30.0 * lplus;
    
    //lx = 0.0014;
    radii[6] = 0.20 * delta;
    radii[7] = 0.20 * delta;
    radii[8] = 0.20 * delta;

    eddy_y_min = (double*)malloc(nregions * sizeof(double));
    eddy_y_max = (double*)malloc(nregions * sizeof(double));
    eddy_y_min[0] = 5.0 * lplus; 
    eddy_y_max[0] = 60.0 * lplus; 
    eddy_y_min[1] = 60.0 * lplus; 
    eddy_y_max[1] = 0.6 * delta;
    eddy_y_min[2] = 0.6 * delta;
    eddy_y_max[2] = y_max + radii[7];

    vols = (double*)malloc(nregions * sizeof(double));
    for (int i{0}; i < nregions; i++){
        vols[i] = std::abs(2 * radii[3*i] * (eddy_y_max[i] - eddy_y_min[i]) * (eddy_z_max - eddy_z_min));
    }

    calc_eddies(eddies, vols, radii, nregions);


    radii[0] = radii[0] / eddy_vel[0];
    radii[3] = radii[3] / eddy_vel[1];
    radii[6] = radii[6] / eddy_vel[2];

    for(int i{0}; i < nregions; i++){
        ops_printf("Region %i -> eddies: %i, vol=%f \n", i, eddies[i], vols[i]);
    }
    // LCG values (change the seed init using LCG), could just set the seed as a random number.
    a = 5;
    c = 3;
    m = pow(2, 29);
    seed_gbl = 2893328493;
    
    // ------------------------------------ Allocate memory for eddy variables -----------------------------

    eddy_pos1 = (double*)malloc(3 * eddies[0] * sizeof(double));
    eddy_r1 = (double*)malloc(3 * eddies[0] * sizeof(double));
    eddy_eps1 = (int*)malloc(3 * eddies[0] * sizeof(int));

    eddy_pos2 = (double*)malloc(3 * eddies[1] * sizeof(double));
    eddy_r2 = (double*)malloc(3 * eddies[1] * sizeof(double));
    eddy_eps2 = (int*)malloc(3 * eddies[1] * sizeof(int));

    eddy_pos3 = (double*)malloc(3 * eddies[2] * sizeof(double));
    eddy_r3 = (double*)malloc(3 * eddies[2] * sizeof(double));
    eddy_eps3 = (int*)malloc(3 * eddies[2] * sizeof(int));

    yinterp = (double*)(y_inp);
    r11interp = (double*)(uu_inp);
    r21interp = (double*)(uv_inp);
    r22interp = (double*)(vv_inp);
    r33interp = (double*)(ww_inp);


    // --------------------------- convert variables to ops const ---------------------------------------------------
    ops_decl_const("u0", 1, "double", &u0);
    ops_decl_const("simulation_time", 1, "double", &simulation_time);
    ops_decl_const("dt", 1, "double", &dt);
    ops_decl_const("delta", 1, "double", &delta);
    ops_decl_const("r_max", 1, "double", &r_max);
    ops_decl_const("ny", 1, "int", &ny);
    ops_decl_const("nz", 1, "int", &nz);
    ops_decl_const("x_min", 1, "double", &x_min);
    ops_decl_const("x_max", 1, "double", &x_max);
    ops_decl_const("x_plane", 1, "double", &x_plane);
    ops_decl_const("y_min", 1, "double", &y_min);
    ops_decl_const("y_max", 1, "double", &y_max);
    ops_decl_const("eddy_y_min", nregions, "double", &eddy_y_min[0]);
    ops_decl_const("eddy_y_max", nregions, "double", &eddy_y_max[0]);
    ops_decl_const("z_min", 1, "double", &z_min);
    ops_decl_const("z_max", 1, "double", &z_max);
    ops_decl_const("z_midpt", 1, "double", &z_midpt);
    ops_decl_const("z_l", 1, "double", & z_l);
    ops_decl_const("eddy_z_min", 1, "double", &eddy_z_min);
    ops_decl_const("eddy_z_max", 1, "double", &eddy_z_max);
    ops_decl_const("Tinf", 1, "double", &Tinf);
    ops_decl_const("Twall", 1, "double", &Twall);
    ops_decl_const("gama", 1, "double", &gama);
    ops_decl_const("Minf", 1, "double", &Minf);
    ops_decl_const("eddies", nregions, "int", &eddies[0]);
    ops_decl_const("radii", 3*nregions, "double", &radii[0]);
    ops_decl_const("eddy_vel", nregions, "double", &eddy_vel[0]);

    
    int eddy_size1[] = {eddies[0], 1};
    int eddy_size2[] = {eddies[1], 1};
    int eddy_size3[] = {eddies[2], 1};
    int eddy_base[] = {0, 0};
    int eddy_pad_minus[] = {0, 0};
    int eddy_pad_plus[] = {0, 0};

    int inlet_size[] = {ny, nz};
    int inlet_base[] = {0, 0};
    int inlet_pad_minus[] = {0, 0};
    int inlet_pad_plus[] = {0, 0};
    

    // inlet plane variables
    double* y_inlet = NULL;
    double* z_inlet = NULL;
    double* uprime = NULL;
    double* vprime = NULL;
    double* wprime = NULL;
    double* Tprime = NULL;
    double* Tbar = NULL;
    double* a11 = NULL;
    double* a21 = NULL;
    double* a22 = NULL;
    double* a31 = NULL;
    double* a32 = NULL;
    double* a33 = NULL;

    // Random number arrays
    // using 1 ops_dat to store all 6 sets of  random numbers: 1x x position, 1x y position, 1x z position, 3x random direction (x, y, z) for each eddy
    int* bulk_rng1 = NULL; // contains 6 random ints per 'grid point'
    int* bulk_rng2 = NULL;
    int* bulk_rng3 = NULL;
    

    ops_block inlet_block = ops_decl_block(2, "inlet");
    ops_block eddy_block = ops_decl_block(2, "eddies");

    // OPS version of eddy variables
    
    double* EDDY_POS1 = NULL;
    double* EDDY_R1 = NULL;
    int* EDDY_EPS1 = NULL;

    double* EDDY_POS2 = NULL;
    double* EDDY_R2 = NULL;
    int* EDDY_EPS2 = NULL;

    double* EDDY_POS3 = NULL;
    double* EDDY_R3 = NULL;
    int* EDDY_EPS3 = NULL;



    ops_dat d_eddy_pos1 = ops_decl_dat(eddy_block, 3, eddy_size1, eddy_base, eddy_pad_minus, eddy_pad_plus, EDDY_POS1, "double", "d_eddy_pos1"); 
    ops_dat d_eddy_r1 = ops_decl_dat(eddy_block, 3, eddy_size1, eddy_base, eddy_pad_minus, eddy_pad_plus, EDDY_R1, "double", "d_eddy_r1");
    ops_dat d_eddy_eps1 = ops_decl_dat(eddy_block, 3, eddy_size1, eddy_base, eddy_pad_minus, eddy_pad_plus, EDDY_EPS1, "int", "d_eddy_eps1");
    ops_dat d_bulk_rng1 = ops_decl_dat(eddy_block, 6, eddy_size1, eddy_base, eddy_pad_minus, eddy_pad_plus, bulk_rng1, "int", "bulk_rng1");

    ops_dat d_eddy_pos2 = ops_decl_dat(eddy_block, 3, eddy_size2, eddy_base, eddy_pad_minus, eddy_pad_plus, EDDY_POS2, "double", "d_eddy_pos2"); 
    ops_dat d_eddy_r2 = ops_decl_dat(eddy_block, 3, eddy_size2, eddy_base, eddy_pad_minus, eddy_pad_plus, EDDY_R2, "double", "d_eddy_r2");
    ops_dat d_eddy_eps2 = ops_decl_dat(eddy_block, 3, eddy_size2, eddy_base, eddy_pad_minus, eddy_pad_plus, EDDY_EPS2, "int", "d_eddy_eps2");
    ops_dat d_bulk_rng2 = ops_decl_dat(eddy_block, 6, eddy_size2, eddy_base, eddy_pad_minus, eddy_pad_plus, bulk_rng2, "int", "bulk_rng2");

    ops_dat d_eddy_pos3 = ops_decl_dat(eddy_block, 3, eddy_size3, eddy_base, eddy_pad_minus, eddy_pad_plus, EDDY_POS3, "double", "d_eddy_pos3"); 
    ops_dat d_eddy_r3 = ops_decl_dat(eddy_block, 3, eddy_size3, eddy_base, eddy_pad_minus, eddy_pad_plus, EDDY_R3, "double", "d_eddy_r3");
    ops_dat d_eddy_eps3 = ops_decl_dat(eddy_block, 3, eddy_size3, eddy_base, eddy_pad_minus, eddy_pad_plus, EDDY_EPS3, "int", "d_eddy_eps3");
    ops_dat d_bulk_rng3 = ops_decl_dat(eddy_block, 6, eddy_size3, eddy_base, eddy_pad_minus, eddy_pad_plus, bulk_rng3, "int", "bulk_rng3");

    // OPS_version of inlet plane variables
    ops_dat d_y_inlet = ops_decl_dat(inlet_block, 1, inlet_size, inlet_base, inlet_pad_minus, inlet_pad_plus, y_inlet, "double", "y_inlet");
    ops_dat d_z_inlet = ops_decl_dat(inlet_block, 1, inlet_size, inlet_base, inlet_pad_minus, inlet_pad_plus, z_inlet, "double", "z_inlet");
    ops_dat d_uprime = ops_decl_dat(inlet_block, 1, inlet_size, inlet_base, inlet_pad_minus, inlet_pad_plus, uprime, "double", "uprime");
    ops_dat d_vprime = ops_decl_dat(inlet_block, 1, inlet_size, inlet_base, inlet_pad_minus, inlet_pad_plus, vprime, "double", "vprime");
    ops_dat d_wprime = ops_decl_dat(inlet_block, 1, inlet_size, inlet_base, inlet_pad_minus, inlet_pad_plus, wprime, "double", "wprime");
    ops_dat d_Tprime = ops_decl_dat(inlet_block, 1, inlet_size, inlet_base, inlet_pad_minus, inlet_pad_plus, Tprime, "double", "Tprime");
    ops_dat d_Tbar = ops_decl_dat(inlet_block, 1, inlet_size, inlet_base, inlet_pad_minus, inlet_pad_plus, Tbar, "double", "Tbar");
    ops_dat d_a11 = ops_decl_dat(inlet_block, 1, inlet_size, inlet_base, inlet_pad_minus, inlet_pad_plus, a11, "double", "a11");
    ops_dat d_a21 = ops_decl_dat(inlet_block, 1, inlet_size, inlet_base, inlet_pad_minus, inlet_pad_plus, a21, "double", "a21");
    ops_dat d_a22 = ops_decl_dat(inlet_block, 1, inlet_size, inlet_base, inlet_pad_minus, inlet_pad_plus, a22, "double", "a22");
    ops_dat d_a31 = ops_decl_dat(inlet_block, 1, inlet_size, inlet_base, inlet_pad_minus, inlet_pad_plus, a31, "double", "a31");
    ops_dat d_a32 = ops_decl_dat(inlet_block, 1, inlet_size, inlet_base, inlet_pad_minus, inlet_pad_plus, a32, "double", "a32");
    ops_dat d_a33 = ops_decl_dat(inlet_block, 1, inlet_size, inlet_base, inlet_pad_minus, inlet_pad_plus, a33, "double", "a33");

    

    // define stencils
    int s1d_00[] = {0};
    ops_stencil S1D_00 = ops_decl_stencil(1, 1, s1d_00, "self1d");
    int eddy_iter_range1[] = {0, eddies[0], 0, 1};
    int eddy_iter_range2[] = {0, eddies[1], 0, 1};
    int eddy_iter_range3[] = {0, eddies[2], 0, 1};

    int s2d_00[] = {0, 0};
    ops_stencil S2D_00 = ops_decl_stencil(2, 1, s2d_00, "self");
    int iter_range[] = {0, ny, 0, nz};

    ops_partition("2D_block_DECOMPSE");

    seed_gbl = (a * seed_gbl + c) % m;
    ops_randomgen_init(seed_gbl, 0);
    ops_fill_random_uniform(d_bulk_rng1);
    ops_fill_random_uniform(d_bulk_rng2);
    ops_fill_random_uniform(d_bulk_rng3);

    // instantiate eddy values
    reg_num = 0;
    
    ops_printf("instantiating eddies\n");
    ops_printf("=====================================\n");
    ops_par_loop(instantiate_eddies, "instantiate_eddies", inlet_block, 2, eddy_iter_range1,
    ops_arg_dat(d_eddy_pos1, 3, S2D_00, "double", OPS_WRITE),
    ops_arg_dat(d_eddy_r1, 3, S2D_00, "double", OPS_WRITE),
    ops_arg_dat(d_eddy_eps1, 3, S2D_00, "int", OPS_WRITE),
    ops_arg_dat(d_bulk_rng1, 6, S2D_00, "int", OPS_READ),
    ops_arg_gbl(&reg_num, 1, "int", OPS_READ));
    reg_num = 1;
    //ops_update_const("reg_num", 1, "int", &reg_num);
    ops_printf("=============================================\n");
    ops_par_loop(instantiate_eddies, "instantiate_eddies", inlet_block, 2, eddy_iter_range2,
    ops_arg_dat(d_eddy_pos2, 3, S2D_00, "double", OPS_WRITE),
    ops_arg_dat(d_eddy_r2, 3, S2D_00, "double", OPS_WRITE),
    ops_arg_dat(d_eddy_eps2, 3, S2D_00, "int", OPS_WRITE),
    ops_arg_dat(d_bulk_rng2, 6, S2D_00, "int", OPS_READ),
    ops_arg_gbl(&reg_num, 1, "int", OPS_READ));

    reg_num = 2;
    //ops_update_const("reg_num", 1, "int", &reg_num);
    ops_printf("=============================================\n");
    ops_par_loop(instantiate_eddies, "instantiate_eddies", inlet_block, 2, eddy_iter_range3,
    ops_arg_dat(d_eddy_pos3, 3, S2D_00, "double", OPS_WRITE),
    ops_arg_dat(d_eddy_r3, 3, S2D_00, "double", OPS_WRITE),
    ops_arg_dat(d_eddy_eps3, 3, S2D_00, "int", OPS_WRITE),
    ops_arg_dat(d_bulk_rng3, 6, S2D_00, "int", OPS_READ),
    ops_arg_gbl(&reg_num, 1, "int", OPS_READ));
    ops_printf("=============================================\n");

    ops_dat_fetch_data(d_eddy_r1, 0, (char*)eddy_r1);
    ops_dat_fetch_data(d_eddy_r2, 0, (char*)eddy_r2);
    ops_dat_fetch_data(d_eddy_r3, 0, (char*)eddy_r3);
    

    ops_par_loop(instantiate_grid, "instantiate_grid", inlet_block, 2, iter_range,
    ops_arg_dat(d_y_inlet, 1, S2D_00, "double", OPS_RW),
    ops_arg_dat(d_z_inlet, 1, S2D_00, "double", OPS_WRITE),
    ops_arg_data(d_Tbar, 1, S2D_00, "double", OPS_WRITE),
    ops_arg_idx());    

    ops_par_loop(instantiate_RST_TBL, "instantiate_RST_TBL", inlet_block, 2, iter_range,
    ops_arg_dat(d_a11, 1, S2D_00, "double", OPS_RW),
    ops_arg_dat(d_a21, 1, S2D_00, "double", OPS_RW),
    ops_arg_dat(d_a22, 1, S2D_00, "double", OPS_RW),
    ops_arg_dat(d_a31, 1, S2D_00, "double", OPS_RW),
    ops_arg_dat(d_a32, 1, S2D_00, "double", OPS_RW),
    ops_arg_dat(d_a33, 1, S2D_00, "double", OPS_RW), 
    ops_arg_dat(d_y_inlet, 1, S2D_00, "double", OPS_READ),
    ops_arg_dat(d_z_inlet, 1, S2D_00, "double", OPS_READ),
    ops_arg_gbl(yinterp, 260, "double", OPS_READ),
    ops_arg_gbl(r11interp, 260, "double", OPS_READ),
    ops_arg_gbl(r21interp, 260, "double", OPS_READ),
    ops_arg_gbl(r22interp, 260, "double", OPS_READ),
    ops_arg_gbl(r33interp, 260, "double", OPS_READ));

    std::string filename;
    
    double ct0, et0;
    double ct1, et1;

    ops_timers(&ct0, &et0);
    ops_printf("%s \n", "------------------------------");

    ops_printf("======================================\n");

    

    for(int i{0}; i < niter; i++){

        simulation_time = simulation_time + dt;
        ops_update_const("simulation_time", 1, "double", &simulation_time);


        ops_fill_random_uniform(d_bulk_rng1);
        ops_fill_random_uniform(d_bulk_rng2);
        ops_fill_random_uniform(d_bulk_rng3);

        reg_num = 0;
        ops_par_loop(convect_eddies, "convect_eddies", eddy_block, 2, eddy_iter_range1,
        ops_arg_dat(d_eddy_pos1, 3, S2D_00, "double", OPS_RW),
        ops_arg_dat(d_eddy_r1, 3, S2D_00, "double", OPS_RW),
        ops_arg_dat(d_eddy_eps1, 3, S2D_00, "int", OPS_WRITE),
        ops_arg_dat(d_bulk_rng1, 6, S2D_00, "int", OPS_READ),
        ops_arg_gbl(&reg_num, 1, "int", OPS_READ));
        
        

        reg_num = 1;
        ops_par_loop(convect_eddies, "convect_eddies", eddy_block, 2, eddy_iter_range2,
        ops_arg_dat(d_eddy_pos2, 3, S2D_00, "double", OPS_RW),
        ops_arg_dat(d_eddy_r2, 3, S2D_00, "double", OPS_RW),
        ops_arg_dat(d_eddy_eps2, 3, S2D_00, "int", OPS_WRITE),
        ops_arg_dat(d_bulk_rng2, 6, S2D_00, "int", OPS_READ),
        ops_arg_gbl(&reg_num, 1, "int", OPS_READ));
        
        

        reg_num = 2;
        ops_par_loop(convect_eddies, "convect_eddies", eddy_block, 2, eddy_iter_range3,
        ops_arg_dat(d_eddy_pos3, 3, S2D_00, "double", OPS_RW),
        ops_arg_dat(d_eddy_r3, 3, S2D_00, "double", OPS_RW),
        ops_arg_dat(d_eddy_eps3, 3, S2D_00, "int", OPS_WRITE),
        ops_arg_dat(d_bulk_rng3, 6, S2D_00, "int", OPS_READ),
        ops_arg_gbl(&reg_num, 1, "int", OPS_READ));
        
        ops_dat_fetch_data(d_eddy_pos1, 0, (char*)eddy_pos1);
        ops_dat_fetch_data(d_eddy_eps1, 0, (char*)eddy_eps1);

        ops_dat_fetch_data(d_eddy_pos2, 0, (char*)eddy_pos2);
        ops_dat_fetch_data(d_eddy_eps2, 0, (char*)eddy_eps2);

        ops_dat_fetch_data(d_eddy_pos3, 0, (char*)eddy_pos3);
        ops_dat_fetch_data(d_eddy_eps3, 0, (char*)eddy_eps3);

        

        ops_par_loop(compute_fluct, "compute_fluct", inlet_block, 2, iter_range,
        ops_arg_dat(d_y_inlet, 1, S2D_00, "double", OPS_READ),
        ops_arg_dat(d_z_inlet, 1, S2D_00, "double", OPS_READ),
        ops_arg_dat(d_a11, 1, S2D_00, "double", OPS_READ),
        ops_arg_dat(d_a21, 1, S2D_00, "double", OPS_READ),
        ops_arg_dat(d_a22, 1, S2D_00, "double", OPS_READ),
        ops_arg_dat(d_a31, 1, S2D_00, "double", OPS_READ),
        ops_arg_dat(d_a32, 1, S2D_00, "double", OPS_READ),
        ops_arg_dat(d_a33, 1, S2D_00, "double", OPS_READ),
        ops_arg_dat(d_uprime, 1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_vprime, 1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_wprime, 1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_Tprime, 1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_Tbar, 1, S2D_00, "double", OPS_READ),
        ops_arg_gbl(eddy_pos1, 3*eddies[0], "double", OPS_READ),
        ops_arg_gbl(eddy_r1, 3*eddies[0], "double", OPS_READ),
        ops_arg_gbl(eddy_eps1, 3*eddies[0], "int", OPS_READ),
        ops_arg_gbl(eddy_pos2, 3*eddies[1], "double", OPS_READ),
        ops_arg_gbl(eddy_r2, 3*eddies[1], "double", OPS_READ),
        ops_arg_gbl(eddy_eps2, 3*eddies[1], "int", OPS_READ),
        ops_arg_gbl(eddy_pos3, 3*eddies[2], "double", OPS_READ),
        ops_arg_gbl(eddy_r3, 3*eddies[2], "double", OPS_READ),
        ops_arg_gbl(eddy_eps3, 3*eddies[2], "int", OPS_READ),
        ops_arg_idx());

        filename = std::string("up-") + std::to_string(i) + ".txt";

        ops_print_dat_to_txtfile(d_uprime, filename.c_str());

    }

    ops_timers(&ct1, &et1);

    ops_printf("time elapsed: %f \n", et1 - et0);

    ops_printf("%s \n", "--------------------");

    ops_exit();
}
