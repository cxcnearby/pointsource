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

TString gtrack_root, geff_root, gdata_root;
TH1F *gsim_nfitc, *gdata_nfitc;
TH1F *ghtrack;
TH3F *ghrecmat;

double Chi2(const double *par);
TH1F *Gethtrack(TString trackroot);
TH3F *Gethrecmat(TString effroot);
void SimNfitcDist(TH1F *ghsim, TH1F *htrack, TH3F *hrecmat,
                  const function<double(double)> &func);
TH1F *DataNfitcDist(TString dataroot);
double CalChi2(TH1F *simnfitcdist, TH1F *datanfitcdist);

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
  gtrack_root = argv[1];
  geff_root = argv[2];
  gdata_root = argv[3];
  ROOT::Math::Minimizer *fit =
      ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  fit->SetMaxFunctionCalls(1000000);
  fit->SetMaxIterations(10000);
  fit->SetTolerance(0.001);
  parinit par;
  gdata_nfitc = DataNfitcDist(gdata_root);
  gsim_nfitc = new TH1F("nfitc", "nfitc dist", kNPmtBin, 0, kPmtRange);
  ROOT::Math::Functor f(&Chi2, par.npar);
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
  SimNfitcDist(gsim_nfitc, ghtrack, ghrecmat, flux);
  return CalChi2(gsim_nfitc, gdata_nfitc);
}

TH1F *Gethtrack(TString trackroot) {
  TFile *ftrack = TFile::Open(trackroot, "READ");
  ghtrack = (TH1F *)ftrack->Get("hzen");
  return ghtrack;
}

TH3F *Gethrecmat(TString effroot) {
  TFile *feff = TFile::Open(effroot, "READ");
  ghrecmat = (TH3F *)feff->Get("recmat0");
  return ghrecmat;
}

void SimNfitcDist(TH1F *ghsim, TH1F *htrack, TH3F *hrecmat,
                  const function<double(double)> &func) {
  int nbinx = hrecmat->GetNbinsX();
  vector<double> Ecuts;
  for (int i = 1; i <= nbinx + 1; ++i) {
    double tmp_cut = hrecmat->GetXaxis()->GetBinLowEdge(i);
    Ecuts.emplace_back(pow(10., tmp_cut));
  }
  vector<double> bin_flux = BinnedIntegratedFlux(IntegratedFlux(func), Ecuts);
  for (int i = 1; i <= kNPmtBin; ++i) {
    double tmp = 0.;
    for (int j = 1; j <= hrecmat->GetNbinsY(); ++j) {
      for (int k = 0; k < nbinx; ++k) {
        tmp += bin_flux[k] * htrack->GetBinContent(j) *
               hrecmat->GetBinContent(k + 1, j, i);
      }
    }
    ghsim->SetBinContent(i, tmp);
  }
}

TH1F *DataNfitcDist(TString dataroot) {
  TFile *fdata = TFile::Open(dataroot, "read");
  TTree *trec = (TTree *)fdata->Get("trec");
  int nfitc;
  trec->SetBranchAddress("nfitc", &nfitc);
  TH1F *simnfitcdist = new TH1F("hnfitc", "nfitc dist", kNPmtBin, 0, kPmtRange);
  Long64_t nentries = trec->GetEntries();
  if (trec->GetBranch("weight")) {
    float weight;
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
