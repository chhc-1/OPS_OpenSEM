void instantiate_grid(ACC<double>& y, ACC<double>& z, const int* idx){
    // pass 1D y, z arrays in, then distribute
    y(0, 0) = y_min + (y_max + r_max - y_min) * (double)(idx[0]) / ny;
    z(0, 0) = z_min + (z_max - z_min) * (double)(idx[1]) / nz;
}

void instantiate_RST(ACC<double>& a11, ACC<double>& a21, ACC<double>& a22, ACC<double>& a31, ACC<double>& a32, ACC<double>& a33){
    a11(0, 0) = u0 * TI;
    a21(0, 0) = 0;
    a22(0, 0) = u0 * TI;
    a31(0, 0) = 0;
    a32(0, 0) = 0;
    a33(0, 0) = u0 * TI;
}

void instantiate_RST_TBL(ACC<double>& a11, ACC<double>& a21, ACC<double>& a22, ACC<double>& a31, ACC<double>& a32, ACC<double>& a33,
                    const ACC<double>& y, const ACC<double>& z, const double* ydata, const double* r11data, const double* r21data, const double* r22data,
                    const double* r33data){
                        
    int idx;
    idx = 0;
    for(int i{1}; i < 260 - 1; i++){
        // assume data passed in is within bounds
        if((y(0, 0) - ydata[i]) < 0){
            idx = i-1;
            break;
        }
    }                   
    double r11temp = (r11data[idx+1] - r11data[idx]) / (ydata[idx+1] - ydata[idx]) * (y(0, 0) - ydata[idx]) + r11data[idx];
    double r21temp = (r21data[idx+1] - r21data[idx]) / (ydata[idx+1] - ydata[idx]) * (y(0, 0) - ydata[idx]) + r21data[idx];
    double r22temp = (r22data[idx+1] - r22data[idx]) / (ydata[idx+1] - ydata[idx]) * (y(0, 0) - ydata[idx]) + r22data[idx];
    double r33temp = (r33data[idx+1] - r33data[idx]) / (ydata[idx+1] - ydata[idx]) * (y(0, 0) - ydata[idx]) + r33data[idx];

    a11(0, 0) = sqrt(r11temp);
    double a11temp = std::abs(a11(0, 0));
    if(a11temp < 0.001){
        a11temp = 0.001;
    }
    a21(0, 0) = r21temp / a11temp;
    a22(0, 0) = sqrt(r22temp - a21(0, 0)*a21(0, 0));
    a31(0, 0) = 0;
    a32(0, 0) = 0;
    a33(0, 0) = sqrt(r33temp - a31(0,0)*a31(0,0) - a32(0,0)*a32(0,0));
}

void convect_eddies(ACC<double>& eddy_pos, ACC<double>& eddy_r, ACC<int>& eddy_eps, const ACC<int>& bulk_rng, const int* reg_num){
    if((simulation_time - eddy_pos(0, 0, 0) - eddy_r(0, 0, 0)) / eddy_r(0, 0, 0) > 1.0){
        eddy_r(0, 0, 0) = radii[3*reg_num[0]];
        eddy_r(1, 0, 0) = radii[3*reg_num[0]+1];
        eddy_r(2, 0, 0) = radii[3*reg_num[0]+2];
        //eddy_pos(0, 0, 0) = simulation_time + ((double)bulk_rng(0,0,0) + 2147483648.0) / (4294967295.0) * eddy_r(0, 0, 0);
        //eddy_pos(1, 0, 0) = eddy_y_min[reg_num[0]] + (eddy_y_max[reg_num[0]] - eddy_y_min[reg_num[0]]) * ((double)bulk_rng(1,0,0) + 2147483648.0) / (4294967295.0);
        //eddy_pos(2, 0, 0) = eddy_z_min + (eddy_z_max - eddy_z_min) * ((double)bulk_rng(2,0,0) + 2147483648.0) / (4294967295.0);
        eddy_pos(0, 0, 0) = simulation_time + ((double)bulk_rng(0,0,0)) / (2147483648.0) * eddy_r(0, 0, 0);
        eddy_pos(1, 0, 0) = eddy_y_min[reg_num[0]] + (eddy_y_max[reg_num[0]] - eddy_y_min[reg_num[0]]) * ((double)bulk_rng(1,0,0)) / 2147483648.0;
        eddy_pos(2, 0, 0) = eddy_z_min + (eddy_z_max - eddy_z_min) * ((double)bulk_rng(2,0,0) ) / (2147483648.0);
        eddy_eps(0, 0, 0) = (bulk_rng(3, 0, 0) < 0) ? -1 : 1;   
        eddy_eps(1, 0, 0) = (bulk_rng(4, 0, 0) < 0) ? -1 : 1;   
        eddy_eps(2, 0, 0) = (bulk_rng(5, 0, 0) < 0) ? -1 : 1;  
        printf("reset eddy in region %i to y=%f, z=%f with bulk rng value 1,2 = %i, %i \n", reg_num[0], eddy_pos(1,0,0), eddy_pos(2,0,0), bulk_rng(1,0,0), bulk_rng(2,0,0));
    }
}

void instantiate_eddies(ACC<double>& eddy_pos, ACC<double>& eddy_r, ACC<int>& eddy_eps, const ACC<int>& bulk_rng, const int* reg_num){
    eddy_r(0, 0, 0) = radii[3*reg_num[0]]; 
    eddy_r(1, 0, 0) = radii[3*reg_num[0]+1];
    eddy_r(2, 0, 0) = radii[3*reg_num[0]+2];
    eddy_pos(0, 0, 0) = simulation_time + (((double)bulk_rng(0, 0, 0) + 2147483648.0) / (4294967295.0) - 1.0) * eddy_r(0, 0, 0);
    eddy_pos(1, 0, 0) = eddy_y_min[reg_num[0]] + ((double)bulk_rng(1, 0, 0) + 2147483648.0) / (4294967295.0) * (eddy_y_max[reg_num[0]] - eddy_y_min[reg_num[0]]);
    // same z max for all eddy regions
    eddy_pos(2, 0, 0) = eddy_z_min + ((double)bulk_rng(2, 0, 0) + 2147483648.0) / (4294967295.0) * (eddy_z_max - eddy_z_min);
    eddy_eps(0, 0, 0) = ((bulk_rng(3, 0, 0) < 0) ? -1 : 1);
    eddy_eps(1, 0, 0) = ((bulk_rng(4, 0, 0) < 0) ? -1 : 1);
    eddy_eps(2, 0, 0) = ((bulk_rng(5, 0, 0) < 0) ? -1 : 1);

}

void compute_fluct(const ACC<double>& y, const ACC<double>& z, const ACC<double>& a11, const ACC<double>& a21, const ACC<double>& a22, const ACC<double>& a31, const ACC<double>& a32, const ACC<double>& a33, 
    ACC<double>& uprime, ACC<double>& vprime, ACC<double>& wprime, const double* eddy_pos1, const double* eddy_r1, const int* eddy_eps1, const double* eddy_pos2, const double* eddy_r2,
    const int* eddy_eps2, const double* eddy_pos3, const double* eddy_r3, const int* eddy_eps3, const int* idx){
    //printf("computing fluct\n");

    double normpos[3] = {0.0, 0.0, 0.0};
    double shape[3] = {0.0, 0.0, 0.0};
    double rsq;
    uprime(0, 0) = 0;
    vprime(0, 0) = 0;
    wprime(0, 0) = 0;

    if (y(0, 0) > 0.0 && y(0, 0) < eddy_y_max[0] + radii[1]){
      for (int i{0}; i < eddies[0]; i++){
          normpos[0] = (simulation_time  - eddy_pos1[3*i] - eddy_r1[3*i]) / (eddy_r1[3*i]);
          normpos[1] = (y(0, 0) - eddy_pos1[3*i+1]) / (eddy_r1[3*i+1]);
          normpos[2] = (z(0, 0) - eddy_pos1[3*i+2]) / (eddy_r1[3*i+2]);
          rsq = normpos[0]*normpos[0] + normpos[1]*normpos[1] + normpos[2]*normpos[2];
  
          if(rsq < 1.0){
              
  
              shape[0] = 1.0/1.7; // scaling factor
              shape[0] *= exp(-0.5*normpos[0]*normpos[0]); // Gaussian
              shape[0] *= exp(-0.5*normpos[1]*normpos[1]); // Gaussian
              shape[0] *= (std::abs(normpos[2]) < 0.001) ? normpos[2] : (1.0 - cos(2*3.14159265*normpos[2])) / (2*3.14158265*normpos[2] * sqrt(0.214)); // 1 - cosine function, change to Gaussian if desired

              shape[1] = 1.0/1.7; // scaling factor
              shape[1] *= exp(-0.5*normpos[0]*normpos[0]); // minus Gaussian due to reversed Cholesky decomp in paper
              shape[1] *= exp(-0.5*normpos[1]*normpos[1]); // Gaussian
              shape[1] *= (std::abs(normpos[2]) < 0.001) ? normpos[2] : (1.0 - cos(2*3.14159265*normpos[2])) / (2*3.14158265*normpos[2] * sqrt(0.214)); // 1 - cosine function
  
              shape[2] = 1.0/1.7; // scaling factor
              shape[2] *= exp(-0.5*normpos[0]*normpos[0]); // Gaussian
              shape[2] *= (std::abs(normpos[1]) < 0.001) ? normpos[1] : (1.0 - cos(2*3.14159265*normpos[1])) / (2*3.14158265*normpos[1] * sqrt(0.214)); // 1 - cosine function
              shape[2] *= exp(-0.5*normpos[2]*normpos[2]); // Gaussian
  
              uprime(0, 0) += a11(0, 0) * eddy_eps1[3*i] * shape[0];
              vprime(0, 0) += (a21(0, 0) * eddy_eps1[3*i] + a22(0, 0) * eddy_eps1[3*i+1]) * shape[1];
              wprime(0, 0) += (a31(0, 0) * eddy_eps1[3*i] + a32(0, 0) * eddy_eps1[3*i+1] + a33(0, 0) * eddy_eps1[3*i+2]) * shape[2];
  
          }
      }
    }
    // hairpin region
    if(y(0, 0) > eddy_y_min[1] - radii[4] && y(0, 0) < eddy_y_max[1] + radii[4]){
      for (int i{0}; i < eddies[1]; i++){
          normpos[0] = (simulation_time - eddy_pos2[3*i] - eddy_r2[3*i]) / (eddy_r2[3*i]);
          normpos[1] = (y(0, 0) - eddy_pos2[3*i+1]) / (eddy_r2[3*i+1]);
          normpos[2] = (z(0, 0) - eddy_pos2[3*i+2]) / (eddy_r2[3*i+2]);
          rsq = normpos[0]*normpos[0] + normpos[1]*normpos[1] + normpos[2]*normpos[2];
  
          
          if(rsq < 1.0){
  
  
              shape[0] = 1.0/1.7; // scaling factor
              shape[0] *= exp(-0.5*normpos[0]*normpos[0]);
              shape[0] *= exp(-0.5*normpos[1]*normpos[1]);
              shape[0] *= (std::abs(normpos[2]) < 0.001) ? normpos[2] : (1.0 - cos(2*3.14159265*normpos[2])) / (2*3.14158265*normpos[2] * sqrt(0.214));
  
              shape[1] = 1.0/1.7; // scaling factor
              shape[1] *= exp(-0.5*normpos[0]*normpos[0]);
              shape[1] *= exp(-0.5*normpos[1]*normpos[1]);
              shape[1] *= (std::abs(normpos[2]) < 0.001) ? normpos[2] : (1.0 - cos(2*3.14159265*normpos[2])) / (2*3.14159265*normpos[2]*sqrt(0.214));
  
              shape[2] = 1.0/1.7; // scaling factor
              shape[2] *= exp(-0.5*normpos[0]*normpos[0]);
              shape[2] *= (std::abs(normpos[1]) < 0.001) ? normpos[1] : (1.0 - cos(2*3.14159265*normpos[1])) / (2*3.14158265*normpos[1] * sqrt(0.214));
              shape[2] *= exp(-0.5*normpos[2]*normpos[2]);
              
              uprime(0, 0) += a11(0, 0) * eddy_eps2[3*i] * shape[0];
              vprime(0, 0) += (a21(0, 0) * eddy_eps2[3*i] + a22(0, 0) * eddy_eps2[3*i+1]) * shape[1];
              wprime(0, 0) += (a31(0, 0) * eddy_eps2[3*i] + a32(0, 0) * eddy_eps2[3*i+1] + a33(0, 0) * eddy_eps2[3*i+2]) * shape[2];
              
          }
      }
    }

    // outer region
    if(y(0, 0) > eddy_y_min[2] - radii[7] && y(0, 0) < eddy_y_max[2] + radii[7]){
      for (int i{0}; i < eddies[2]; i++){
          normpos[0] = (simulation_time - eddy_pos3[3*i] - eddy_r3[3*i]) / (eddy_r3[3*i]);
          normpos[1] = (y(0, 0) - eddy_pos3[3*i+1]) / (eddy_r3[3*i+1]);
          normpos[2] = (z(0, 0) - eddy_pos3[3*i+2]) / (eddy_r3[3*i+2]);
          rsq = normpos[0]*normpos[0] + normpos[1]*normpos[1] + normpos[2]*normpos[2];
  
  
          
          if(rsq < 1.0){
              shape[0] = 1.0/1.4; // scaling factor
              shape[0] *= exp(-0.5*normpos[0]*normpos[0]);
              shape[0] *= exp(-0.5*normpos[1]*normpos[1]);
              shape[0] *= exp(-0.5*normpos[2]*normpos[2]);
  
              shape[1] = 1.0/1.4; // scaling factor
              shape[1] *= exp(-0.5*normpos[0]*normpos[0]);
              shape[1] *= exp(-0.5*normpos[1]*normpos[1]);
              shape[1] *= exp(-0.5*normpos[2]*normpos[2]);
  
              shape[2] = 1.0/1.4; // scaling factor
              shape[2] *= exp(-0.5*normpos[0]*normpos[0]);
              shape[2] *= exp(-0.5*normpos[1]*normpos[1]);
              shape[2] *= exp(-0.5*normpos[2]*normpos[2]);
  
              uprime(0, 0) += a11(0, 0) * eddy_eps3[3*i] * shape[0];
              vprime(0, 0) += (a21(0, 0) * eddy_eps3[3*i] + a22(0, 0) * eddy_eps3[3*i+1]) * shape[1];
              wprime(0, 0) += (a31(0, 0) * eddy_eps3[3*i] + a32(0, 0) * eddy_eps3[3*i+1] + a33(0, 0) * eddy_eps3[3*i+2]) * shape[2];

          }
      }
    }
}

