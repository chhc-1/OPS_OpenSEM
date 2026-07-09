
void calc_eddies(int& eddies, const double& vol, const double& rep_radius){
    eddies = trunc(vol / pow(rep_radius, 3));
}