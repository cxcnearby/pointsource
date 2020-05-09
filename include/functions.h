#ifndef POINTSOURCE_INCLUDE_FUNCTIONS_H_
#define POINTSOURCE_INCLUDE_FUNCTIONS_H_

#include <vector>

double flux(const double E, const int type);
double integrated_flux(const double E1, const double E2);
std::vector<double> binned_integrated_flux(const int type);
std::vector<double> duration_of_zenith_bin();
std::vector<double> stripratio(const int type, const double zen_bin_width,
                               const int n_zen_bin);

#endif
