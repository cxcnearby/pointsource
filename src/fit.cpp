#include "constant.h"
#include "functions.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TFile.h"
#include "TH1.h"
#include "TH3.h"
#include "TTree.h"

using namespace std;

TString track_root, eff_root, data_root;

double Chi2(const double *par);
TH1F *SimNfitcDist(TString trackroot, TString effroot,
                   const function<double(double)> &func);
TH1F *DataNfitcDist(TString dataroot);
double CalChi2(TH1F *simnfitcdist, TH1F *datanfitcdist);
double test(const double *par) { return par[0] + par[1] + par[2] + par[3]; }

struct parinit {
  int npar = 4;
  double var[4] = {1., 1., -4., -1.};
  double step[4] = {0.01, 0.01, 0.01, 0.01};
  string parname[4] = {"F0", "E0", "A", "B"};
};

int main(int argc, char *argv[]) {
  if (argc < 4) {
    cout << argv[0] << "  trackroot  effroot  dataroot" << endl;
    exit(0);
  }
  track_root = argv[1];
  eff_root = argv[2];
  data_root = argv[3];
  // ROOT::Minuit2::Minuit2Minimizer *fit =
  //     new ROOT::Minuit2::Minuit2Minimizer(ROOT::Minuit2::kMigrad);
  ROOT::Math::Minimizer *fit =
      ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  printf("test1\n");
  fit->SetMaxFunctionCalls(1000000);
  fit->SetMaxIterations(10000);
  fit->SetTolerance(0.001);
  parinit par;
  ROOT::Math::Functor f(&Chi2, par.npar);
  printf("test2\n");
  fit->SetFunction(f);
  for (int i = 0; i < par.npar; ++i) {
    fit->SetVariable(i, par.parname[i], par.var[i], par.step[i]);
  }
  fit->Minimize();
  return 0;
}

double Chi2(const double *par) {
  const double F0 = par[0] * 1.e-9;
  const double E0 = par[1];
  const double A = par[2];
  const double B = par[3];
  auto flux = [&](const double E) {
    return F0 * pow(E / E0, A + B * log(E / E0));
  };
  TH1F *sim_nfitc = SimNfitcDist(track_root, eff_root, flux);
  TH1F *data_nfitc = DataNfitcDist(data_root);
  return CalChi2(sim_nfitc, data_nfitc);
}

TH1F *SimNfitcDist(TString trackroot, TString effroot,
                   const function<double(double)> &func) {
  vector<double> point_time_zen = PointDurationOfZenithBin(trackroot);
  TFile *fsim = TFile::Open(effroot, "READ");
  TH3F *h_e0_zenmc_nfitc = (TH3F *)fsim->Get("h_e0_zenmc_nfitc_0");
  int nbinx = h_e0_zenmc_nfitc->GetNbinsX();
  vector<double> Ecuts;
  for (int i = 1; i <= nbinx + 1; ++i) {
    double tmp_cut = h_e0_zenmc_nfitc->GetXaxis()->GetBinLowEdge(i);
    Ecuts.emplace_back(pow(10., tmp_cut));
  }
  vector<double> bin_flux = BinnedIntegratedFlux(IntegratedFlux(func), Ecuts);
  TH1F *simnfitcdist = new TH1F("nfitc", "nfitc dist", kNPmtBin, 0, kPmtRange);
  for (int i = 1; i <= kNPmtBin; ++i) {
    double tmp = 0.;
    for (int j = 0; j < kZenRange; ++j) {
      for (int k = 0; k < nbinx; ++k) {
        tmp += bin_flux[k] * point_time_zen[j] *
               h_e0_zenmc_nfitc->GetBinContent(k + 1, j + 1, i);
      }
    }
    simnfitcdist->SetBinContent(i, tmp);
  }
  return simnfitcdist;
}

TH1F *DataNfitcDist(TString dataroot) {
  TFile *fdata = TFile::Open(dataroot, "read");
  TTree *trec = (TTree *)fdata->Get("trec");
  int nfitc;
  trec->SetBranchAddress("nfitc", &nfitc);
  TH1F *simnfitcdist = new TH1F("nfitc", "nfitc dist", kNPmtBin, 0, kPmtRange);
  Long64_t nentries = trec->GetEntries();
  if (trec->GetBranch("weight")) {
    double weight;
    trec->SetBranchAddress("weight", &weight);
    for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
      trec->GetEntry(ientry);
      simnfitcdist->Fill(nfitc, weight);
    }
  } else {
    for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
      trec->GetEntry(ientry);
      simnfitcdist->Fill(nfitc);
    }
  }
  return simnfitcdist;
}

double CalChi2(TH1F *simnfitcdist, TH1F *datanfitcdist) {
  double chi2 = 0.;
  for (int i = 1; i <= simnfitcdist->GetNbinsX(); ++i) {
    double yd = datanfitcdist->GetBinContent(i);
    double ys = simnfitcdist->GetBinContent(i);
    double err = simnfitcdist->GetBinError(i);
    if (err < 1.e-8)
      continue;
    chi2 += pow((yd - ys) / err, 2);
  }
  return chi2;
}
