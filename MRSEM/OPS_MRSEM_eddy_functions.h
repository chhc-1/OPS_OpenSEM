
void calc_eddies(int* eddies, const double* vol, const double* radii, const int& nregions){
    for(int i{0}; i < nregions; i++){
        eddies[i] = trunc(vol[i]/(radii[3*i]*radii[3*i+1]*radii[3*i+2]));
    }
}