#ifndef POINTSOURCE_INCLUDE_FUNCTIONS_H_
#define POINTSOURCE_INCLUDE_FUNCTIONS_H_

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "TTree.h"

#include <functional>
#include <vector>

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
template <typename T>
TH1F *FillHist1D(const TString &rootname, const TString &treename,
                 const TString &histname, const TString &histtitle,
                 const TString &xname, T &x, const int &xbins,
                 const double &xlow, const double &xup) {
  TH1F *h1 = new TH1F(histname, histtitle, xbins, xlow, xup);
  TFile *f = TFile::Open(rootname, "read");
  TTree *t = (TTree *)f->Get(treename);
  t->SetBranchAddress(xname, &x);
  Long64_t nentries = t->GetEntries();
  for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
    t->GetEntry(ientry);
    h1->Fill(log10(x / 1000.));
  }
  f->Close();
  return h1;
}
template <typename T1, typename T2>
TH2F *FillHist2D(const TString &rootname, const TString &treename,
                 const TString &histname, const TString &histtitle,
                 const TString &xname, T1 &x, const int &xbins,
                 const double &xlow, const double &xup, const TString &yname,
                 T2 &y, const int &ybins, const double &ylow,
                 const double &yup) {
  TH2F *h2 = new TH2F(histname, histtitle, xbins, xlow, xup, ybins, ylow, yup);
  TFile *f = TFile::Open(rootname, "read");
  TTree *t = (TTree *)f->Get(treename);
  t->SetBranchAddress(xname, &x);
  t->SetBranchAddress(yname, &y);
  Long64_t nentries = t->GetEntries();
  for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
    t->GetEntry(ientry);
    h2->Fill(log10(x / 1000.), y);
  }
  f->Close();
  return h2;
}
template <typename T1, typename T2, typename T3>
TH3F *FillHist3D(const TString &rootname, const TString &treename,
                 const TString &histname, const TString &histtitle,
                 const TString &xname, T1 &x, const int &xbins,
                 const double &xlow, const double &xup, const TString &yname,
                 T2 &y, const int &ybins, const double &ylow, const double &yup,
                 const TString &zname, T3 &z, const int &zbins,
                 const double &zlow, const double &zup) {
  TH3F *h3 = new TH3F(histname, histtitle, xbins, xlow, xup, ybins, ylow, yup,
                      zbins, zlow, zup);
  TFile *f = TFile::Open(rootname, "read");
  TTree *t = (TTree *)f->Get(treename);
  t->SetBranchAddress(xname, &x);
  t->SetBranchAddress(yname, &y);
  t->SetBranchAddress(zname, &z);
  Long64_t nentries = t->GetEntries();
  for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
    t->GetEntry(ientry);
    h3->Fill(log10(x / 1000.), y, z);
  }
  f->Close();
  return h3;
}

#endif
