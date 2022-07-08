auto id_pair(RVecF v, RVecI idx);
auto id_filter(RVecF v, RVecI idx);
auto vector_division(RVecF num, RVecF denom);
auto flavour_match(RVecF arr, RVecI flavs, int type);
auto flavour_match(RVecF arr, RVecI flavs, RVecI type);
auto dR_filt(RVecF eta1, RVecF phi1, RVecF eta2, RVecF phi2);
auto leading_two(RVecF v);
auto leading_twoarg(RVecF arr);
double Median(const TH1D * h1);
RVecD fill_median(TH2 *h, int nxbins); 
RVecD fill_gmean(TH2 *h, int nxbins, double Nsigma);
auto veto_map_filt(TH2 *h, RVecF eta, RVecF rho, RVecF input);