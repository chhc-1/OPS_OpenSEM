

void instantiate_grid(ACC<double>& y, ACC<double>& z, const int* idx){
    // pass 1D y, z arrays in, then distribute
    y(0, 0) = y_min + (y_max + r_max - y_min) * (double)(idx[0]) / ny;
    z(0, 0) = z_min - r_max + (z_max + r_max - z_min + r_max) * (double)(idx[1]) / nz;
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

void convect_eddies(ACC<double>& x, ACC<double>& y, ACC<double>& z, ACC<double>& radius, ACC<double>& increment,
                    ACC<int>& eps_x, ACC<int>& eps_y, ACC<int>& eps_z, const ACC<int>& y_rng, const ACC<int>& z_rng,
                    const ACC<int>& eps_x_rng, const ACC<int>& eps_y_rng, const ACC<int>& eps_z_rng){
    x(0, 0) += increment(0, 0);

    
    double xnorm = 0.0;
    double shift = 0.0;
    double shifted_yrng = 0.0;
    
    if(x(0, 0) > x_max){
        xnorm = (double)y_rng(0, 0)/1073741824.0;
        shift = (xnorm + k) / b;
        shifted_yrng = rng_minval/pow(shift, powfact);

        x(0, 0) -= 2 * r_max;
        y(0, 0) = eddy_y_min + (eddy_y_max - eddy_y_min) * (shifted_yrng);
        //y(0, 0) = eddy_y_min + (eddy_y_max - eddy_y_min) * (((double)y_rng(0, 0) + 2147483648.0) / (4294967295.0));
        z(0, 0) = eddy_z_min + (eddy_z_max - eddy_z_min) * (((double)z_rng(0, 0)) / (2147483648.0));
        eps_x(0, 0) = ((eps_x_rng(0, 0) < 1073741824) ? -1 : 1);
        eps_y(0, 0) = ((eps_y_rng(0, 0) < 1073741824) ? -1 : 1);
        eps_z(0, 0) = ((eps_z_rng(0, 0) < 1073741824) ? -1 : 1);
        radius(0, 0) = (y(0, 0) < dmin) ? dmin : 0.41 * y(0,0);
        increment(0, 0) = (y(0, 0) < delta * 0.01) ? pow(0.01, 0.16666667) * dt : (pow(y(0, 0)/delta, 0.16666667) * dt);
        //printf("reset an eddy. New y pos: %f. New radius: %f New increment: %f. Shifted rng value: %f. Shift value: %f. xnorm value: %f. yrng value: %i \n", y(0, 0), radius(0, 0), increment(0, 0), shifted_yrng, shift, xnorm, y_rng(0, 0));
    }
}

void instantiate_eddies(ACC<double>& x, ACC<double>& y, ACC<double>& z, ACC<double>& radius, ACC<double>& increment,
    ACC<int>& eps_x, ACC<int>& eps_y, ACC<int>& eps_z, const ACC<int>& x_rng, const ACC<int>& y_rng, const ACC<int>& z_rng, 
    const ACC<int>& eps_x_rng, const ACC<int>& eps_y_rng, const ACC<int>& eps_z_rng){
    
    double xnorm = 0.0;
    double shift = 0.0;
    double shifted_yrng = 0.0;


    xnorm = (double)y_rng(0, 0)/1073741824.0;
    shift = (xnorm + k) / b;
    shifted_yrng = rng_minval/pow(shift, powfact);
    

    x(0, 0) = x_min + ((double)x_rng(0, 0)) / (2147483648.0) * (x_max - x_min);
    y(0, 0) = eddy_y_min + (eddy_y_max - eddy_y_min) * (shifted_yrng);
    z(0, 0) = eddy_z_min + ((double)z_rng(0, 0)) / (2147483648.0) * (eddy_z_max - eddy_z_min);
    radius(0, 0) = (y(0, 0) < dmin) ? dmin :  y(0, 0) * 0.41;
    increment(0, 0) = (y(0, 0) < delta * 0.01) ? pow(0.01, 0.16666667) * dt : (pow(y(0, 0)/delta, 0.16666667) * dt);
    eps_x(0, 0) = ((eps_x_rng(0, 0) < 1073741824) ? -1 : 1);
    eps_y(0, 0) = ((eps_y_rng(0, 0) < 1073741824) ? -1 : 1);
    eps_z(0, 0) = ((eps_z_rng(0, 0) < 1073741824) ? -1 : 1);
}

void instantiate_Tbar(const ACC<double>& y, ACC<double>& Tbar){
    double u = 0.0;
    u = (y(0,0) < delta) ? u0 * pow(y(0,0)/delta, 0.166667): u0;

    Tbar(0, 0) = Twall + (Tinf-Twall)*u*u;
}

void compute_fluct(const ACC<double>& y, const ACC<double>& z, const ACC<double>& a11, const ACC<double>& a21, const ACC<double>& a22, const ACC<double>& a31, const ACC<double>& a32, const ACC<double>& a33, 
    ACC<double>& uprime, ACC<double>& vprime, ACC<double>& wprime, ACC<double>& Tprime, const ACC<double>& Tbar, const double* x_gbl, const double* y_gbl, const double* z_gbl,
    const double* r_gbl, const int* eps_x_gbl, const int* eps_y_gbl, const int* eps_z_gbl){
    
    double shape = 0.0;
    double rsq;
    //double shape_temp = 0.0;
    //int idxA;
    uprime(0, 0) = 0;
    vprime(0, 0) = 0;
    wprime(0, 0) = 0;

    for (int i{0}; i < eddies; i++){
        rsq = x_gbl[i]*x_gbl[i] + (y_gbl[i]-y(0,0))*(y_gbl[i]-y(0,0)) + (z_gbl[i]-z(0,0))*(z_gbl[i]-z(0,0));
        if(rsq < r_gbl[i]*r_gbl[i]){
            shape = (std::abs(x_gbl[i] - x_plane) < r_gbl[i]) ? exp(-0.5*(x_gbl[i]-x_plane)*(x_gbl[i]-x_plane)/(r_gbl[i]*r_gbl[i])) : 0.0;
            shape *= 1/1.5829045; //1/2.2385651; //1/5.01117;
            shape *= exp(-0.5*(y_gbl[i]-y(0,0))*(y_gbl[i]-y(0,0))/(r_gbl[i]*r_gbl[i]));
            shape *= exp(-0.5*(z_gbl[i]-z(0,0))*(z_gbl[i]-z(0,0))/(r_gbl[i]*r_gbl[i]));
            

            uprime(0, 0) += a11(0, 0) * eps_x_gbl[i] * shape;
            vprime(0, 0) += a21(0, 0) * eps_x_gbl[i] * shape + a22(0, 0) * eps_y_gbl[i] * shape;
            wprime(0, 0) += a31(0, 0) * eps_x_gbl[i] * shape + a32(0, 0) * eps_y_gbl[i] * shape + a33(0, 0) * eps_z_gbl[i] * shape;
        }
    }
    uprime(0, 0) = uprime(0, 0) * vsf;
    vprime(0, 0) = vprime(0, 0) * vsf;
    wprime(0, 0) = wprime(0, 0) * vsf;
    Tprime(0, 0) = -Tbar(0,0)*(gama-1.0)*Minf*Minf*uprime(0,0)/u0;
}

