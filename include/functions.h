#ifndef POINTSOURCE_INCLUDE_FUNCTIONS_H_
#define POINTSOURCE_INCLUDE_FUNCTIONS_H_

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "TTree.h"

#include <functional>
#include <regex>
#include <vector>

int GetTypeFromDir(const std::string &dir);
int TypeCor2TypeId(const int typecor);
double GammaFlux(const double E);
double CRFlux(const int type, const double E);
std::function<double(double)> Flux(const int type);
std::function<double(double, double)>
IntegratedFlux(const std::function<double(double)> &func);
std::vector<double>
BinnedIntegratedFlux(const std::function<double(double, double)> &func,
                     std::vector<double> energy_bin);
std::vector<double> StripArea(const double z1, const double z2,
                              const double zwidth);
std::vector<double> PointDurationOfZenithBin(const TString trackstatroot);
std::vector<double> InWindowDurationOfZenithBin(const TString trackstatroot,
                                                const double window_radius,
                                                const double direction_error);

#endif
