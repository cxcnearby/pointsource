#include "functions.h"
#include "constant.h"

#include <iostream>
#include <math.h>

int TypeCor2TypeId(const int typecor) {
  int i_type = -999;
  switch (typecor) {
  case 0:
    i_type = 0;
    break;
  case 14:
    i_type = 1;
    break;
  case 402:
    i_type = 2;
    break;
  case 1407:
    i_type = 3;
    break;
  case 2513:
    i_type = 4;
    break;
  case 5626:
    i_type = 5;
    break;
  default:
    break;
  }
  return i_type;
}

double GammaFlux(const double E) {
  double F0 = 2.35e-9; // TeV m^-2 s^-1
  double E0 = 7.;      // TeV
  double A = -2.79;
  double B = -0.10;
  return F0 * pow(E / E0, A + B * log(E / E0));
}

double CRFlux(const int type, const double E) {
  double kCRpar[][4][3] = {{{7860., 1.66, 4.e3},
                            {20., 1.4, 3.e4},
                            {1.7, 1.4, 2.e6},
                            {200., 1.6, 6.e7}},
                           {{3550., 1.58, 4.e3},
                            {20., 1.4, 3.e4},
                            {1.7, 1.4, 2.e6},
                            {0.0, 0., 6.e7}},
                           {{2200., 1.63, 4.e3},
                            {13.4, 1.4, 3.e4},
                            {1.14, 1.4, 2.e6},
                            {0., 0., 6.e7}},
                           {{1430., 1.67, 4.e3},
                            {13.4, 1.4, 3.e4},
                            {1.14, 1.4, 2.e6},
                            {0., 0., 6.e7}},
                           {{2120., 1.63, 4.e3},
                            {13.4, 1.4, 3.e4},
                            {1.14, 1.4, 2.e6},
                            {0., 0., 6.e7}}};
  double z[] = {1., 2., 7., 13., 26.};
  auto flux = [=](const int i_type) {
    double flux_part = 0.;
    for (unsigned int i = 0; i < sizeof(kCRpar[0][0]) / sizeof(kCRpar[0][0][0]);
         ++i) {
      flux_part += 1.e3 * kCRpar[i_type][i][0] *
                   pow(1000. * E, -1. - kCRpar[i_type][i][1]) *
                   exp(-E / (z[i_type] * kCRpar[i_type][i][2]));
    }
    return flux_part;
  };
  return flux(TypeCor2TypeId(type));
}

std::function<double(double)> Flux(const int type) {
  if (type == 0) {
    return [=](double E) { return GammaFlux(E); };
  } else {
    return [=](double E) { return CRFlux(type, E); };
  }
}

std::function<double(double, double)>
IntegratedFlux(const std::function<double(double)> &func) {
  auto integratedflux = [=](double E1, double E2) {
    double N = 0;
    int n = 500;
    double dindex = pow(E2 / E1, 1. / n);
    double E = E1;
    for (int i = 0; i < n; ++i) {
      double dE = E * (dindex - 1.);
      double dN = (func(E) + 4. * func(E + dE / 2.) + func(E + dE)) * dE / 6.;
      N += dN;
      E += dE;
    }
    return N;
  };
  return integratedflux;
}

std::vector<double>
BinnedIntegratedFlux(const std::function<double(double, double)> &func,
                     std::vector<double> energy_bin) {
  // unequal-bin-width compound Simpson method.
  std::vector<double> binned_integrated_flux;
  for (unsigned int i = 0; i < energy_bin.size() - 1; ++i) {
    double N = func(energy_bin[i], energy_bin[i + 1]);
    binned_integrated_flux.emplace_back(N);
  }
  return binned_integrated_flux;
}

std::vector<double> StripArea(const double z1, const double z2,
                              const double zwidth) {
  std::vector<double> striparea;
  int n_z = (z2 - z1) / zwidth;
  for (int i = 0; i < n_z; ++i) {
    double tmp_area =
        2. * PI *
        (cos((z1 + i * zwidth) * D2R) - cos((z1 + (i + 1) * zwidth) * D2R));
    striparea.emplace_back(tmp_area);
  }
  return striparea;
}

std::vector<double> PointDurationOfZenithBin(const TString trackstatroot) {
  TFile *fcrab = TFile::Open(trackstatroot, "read");
  TH1F *hzen = (TH1F *)fcrab->Get("hzen");
  const int kNzen = hzen->GetNbinsX();
  std::vector<double> duration;
  for (int i = 1; i <= kNzen; ++i) {
    duration.emplace_back(hzen->GetBinContent(i));
  }
  fcrab->Close();
  return duration;
}

std::vector<double> InWindowDurationOfZenithBin(const TString trackstatroot,
                                                const double window_radius,
                                                const double direction_error) {
  const double a = (window_radius + direction_error) * D2R;
  std::vector<double> duration;
  std::vector<double> point_duration = PointDurationOfZenithBin(trackstatroot);
  const int kNzen = point_duration.size();
  for (int i = 0; i < kNzen; ++i) {
    double eff_time = 0.;
    const double b = (i + 0.5) * kZenBinWidth * D2R;
    for (int j = 0; j < kNzen; ++j) {
      const double c = (j + 0.5) * kZenBinWidth * D2R;
      double A;
      if (fabs(c - b) > a) // two circles have no overlap.
        continue;
      if (c + b < a) { // window circle contains strip circle.
        A = PI;
      } else { // two circles have some overlap.
        double cccss = (cos(a) - cos(b) * cos(c)) / (sin(b) * sin(c));
        if (cccss > 1.) // to ensure fabs(cccss) < 1
          cccss = 1.;
        if (cccss < -1.)
          cccss = -1.;
        A = acos(cccss);
      }
      double segment_ratio = A / PI;
      double time_of_crab_inbin = point_duration[j];
      eff_time += segment_ratio * time_of_crab_inbin;
    }
    duration.emplace_back(eff_time);
  }
  return duration;
} // cosa = cosb cosc + sinb sinc cosA, derive A!
