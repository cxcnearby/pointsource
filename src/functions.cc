#include "functions.h"
#include "constant.h"

#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include <math.h>

double flux(const double E, const int type) {
  double flux = 0.;
  switch (type) {
  case 0: {
    // Crab spectrum
    double F0 = 2.35e-9; // TeV m^-2 s^-1
    double E0 = 7.;      // TeV
    double A = -2.79;
    double B = -0.10;
    flux = F0 * pow(E / E0, A + B * log(E / E0));
    break;
  }
  case 14: {
    // Gaisser model
    double par[4][3] = {{7860., 1.66, 4.e3},
                        {20., 1.4, 3.e4},
                        {1.7, 1.4, 2.e6},
                        {200., 1.6, 6.e7}};
    double z = 1.;
    for (int i = 0; i < 3; ++i) {
      flux += 1.e3 * par[i][0] * pow(1000. * E, -1. - par[i][1]) *
              exp(-E / (z * par[i][2])); // TeV m^-2 s^-1 sr^-1
    }
    break;
  }
  case 402: {
    double par[4][3] = {{3550., 1.58, 4.e3},
                        {20., 1.4, 3.e4},
                        {1.7, 1.4, 2.e6},
                        {0., 0., 6.e7}};
    double z = 2.;
    for (int i = 0; i < 3; ++i) {
      flux += 1.e3 * par[i][0] * pow(1000. * E, -1. - par[i][1]) *
              exp(-E / (z * par[i][2]));
    }
    break;
  }
  case 1407: {
    double par[4][3] = {{2200., 1.63, 4.e3},
                        {13.4, 1.4, 3.e4},
                        {1.14, 1.4, 2.e6},
                        {0., 0., 6.e7}};
    double z = 7.;
    for (int i = 0; i < 3; ++i) {
      flux += 1.e3 * par[i][0] * pow(1000. * E, -1. - par[i][1]) *
              exp(-E / (z * par[i][2]));
      // cout << E << "  "<< par[i][0] << "  "<< pow(1000. * E, -1. - par[i][1])
      // << "  " << exp(-E / (z * par[i][2])) <<"  "<<flux << endl;
    }
    break;
  }
  case 2513: {
    double par[4][3] = {{1430., 1.67, 4.e3},
                        {13.4, 1.4, 3.e4},
                        {1.14, 1.4, 2.e6},
                        {0., 0., 6.e7}};
    double z = 13.;
    for (int i = 0; i < 3; ++i) {
      flux += 1.e3 * par[i][0] * pow(1000. * E, -1. - par[i][1]) *
              exp(-E / (z * par[i][2]));
    }
    break;
  }
  case 5626: {
    double par[4][3] = {{2120., 1.63, 4.e3},
                        {13.4, 1.4, 3.e4},
                        {1.14, 1.4, 2.e6},
                        {0., 0., 6.e7}};
    double z = 26.;
    for (int i = 0; i < 3; ++i) {
      flux += 1.e3 * par[i][0] * pow(1000. * E, -1. - par[i][1]) *
              exp(-E / (z * par[i][2]));
    }
    break;
  }
  default:
    std::cerr << "Wrong Particle Type" << std::endl;
    exit(1);
  }
  return flux;
}

double integrated_flux(const int type, const double E1, const double E2) {
  double N = 0;
  int n = 500;
  double dindex = pow(E2 / E1, 1. / n);
  double E = E1;
  for (int i = 0; i < n; ++i) {
    double dE = E * (dindex - 1.);
    double dN =
        (flux(E, type) + 4. * flux(E + dE / 2., type) + flux(E + dE, type)) *
        dE / 6.;
    N += dN;
    E += dE;
  }
  return N;
}

std::vector<double> binned_integrated_flux(const int type) {
  // unequal-bin-width compound Simpson method.
  std::vector<double> binned_integrated_flux;
  double energy_bin[] = {0.01, 0.1, 1., 1.e1, 1.e2, 1.e3};
  int bin_number = sizeof(energy_bin) / sizeof(energy_bin[0]) - 1;
  for (int i = 0; i < bin_number; ++i) {
    double N = integrated_flux(type, energy_bin[i], energy_bin[i + 1]);
    binned_integrated_flux.emplace_back(N);
  }
  return binned_integrated_flux;
}

std::vector<double> duration_of_zenith_bin() {
  TFile *fcrab = TFile::Open("crab_zen_dist.root", "read");
  TH1F *hzen = (TH1F *)fcrab->Get("hzen");
  const int kNzen = hzen->GetNbinsX();
  std::vector<double> duration;
  for (int i = 0; i < kNzen; ++i) {
    duration.emplace_back(hzen->GetBinContent(i));
  }
  fcrab->Close();
  return duration;
}

std::vector<double> stripratio(const int type, const double zen_bin_width,
                               const int n_zen_bin) {
  std::vector<double> strip;
  for (int i = 0; i < n_zen_bin; ++i) {
    double totalarea = 2. * PI * (1. - cos(60. * D2R));
    double tmp_stripratio;
    if (type != 0) {
      tmp_stripratio =
          1. / (totalarea * cos((i + 0.5) * zen_bin_width *
                                D2R)); // 1/(A_tot*cos(theta)) for CR
    } else {
      double striparea =
          2. * PI *
          (cos(i * zen_bin_width * D2R) - cos((i + 1) * zen_bin_width * D2R));
      tmp_stripratio =
          striparea / (totalarea * cos((i + 0.5) * zen_bin_width * D2R));
    } // A_strip/(A_tot*cos(theta)) for Crab;
    strip.emplace_back(tmp_stripratio);
  }
  return strip;
}
