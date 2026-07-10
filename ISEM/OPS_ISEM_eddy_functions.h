void calc_eddies(int& eddies, const double& vol, const double& rep_radius){
    eddies = trunc(vol / (rep_radius*rep_radius*rep_radius));
}