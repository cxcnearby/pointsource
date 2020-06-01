#include "constant.h"
#include "functions.h"

#include <iostream>
#include <memory>

using namespace std;

void CalcEffareaResponse(TString inputroot);
inline double DivideNoNan(double up, double down) {
  return ((down < 1.e-8) ? 0 : up / down);
};

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << argv[0] << "  weight_added_root_file" << endl;
    exit(0);
  }
  TString input = argv[1];
  CalcEffareaResponse(input);
  return 0;
}

void CalcEffareaResponse(TString inputroot) {
  shared_ptr<TFile> fsim(TFile::Open(inputroot, "read"));
  shared_ptr<TTree> trec((TTree *)fsim->Get("trec"));
  float e0, zenmc, compactness, weight;
  int nfitc, type;
  trec->SetBranchAddress("e0", &e0);
  trec->SetBranchAddress("zenmc", &zenmc);
  trec->SetBranchAddress("compactness", &compactness);
  trec->SetBranchAddress("weight", &weight);
  trec->SetBranchAddress("nfitc", &nfitc);
  trec->SetBranchAddress("type", &type);
  TString outroot = inputroot.ReplaceAll(".root", "_effdir.root");
  shared_ptr<TFile> feff(TFile::Open(outroot, "RECREATE"));
  vector<TH2F> vhrec_eff, vhsel_eff;
  vector<TH3F> vhrec_mat, vhsel_mat;
  vector<shared_ptr<TH2F>> vh2d;
  for (auto &&v : kTypeCode) {
    TH2F hrec_eff(Form("receff%d", v), Form("receff%d--energy x zenith", v),
                  kNEnergyBin, -2, 3, kNZenBin, 0,
                  kZenRange); // effarea of reconstructed events
    TH2F hsel_eff(Form("seleff%d", v), Form("seleff%d--energy x zenith", v),
                  kNEnergyBin, -2, 3, kNZenBin, 0,
                  kZenRange); // effarea of events after selection
    TH3F hrec_mat(Form("recmat%d", v),
                  Form("recmat%d--energy x zenmc x nfitc", v), kNEnergyBin, -2,
                  3, kNZenBin, 0, kZenRange, kNPmtBin, 0,
                  kPmtRange); // response of reconstructed events
    TH3F hsel_mat(Form("selmat%d", v),
                  Form("selmat%d--energy x zenmc x nfitc", v), kNEnergyBin, -2,
                  3, kNZenBin, 0, kZenRange, kNPmtBin, 0,
                  kPmtRange); // response of events after selection
    shared_ptr<TH2F> h2d((TH2F *)fsim->Get(Form("h_e0_zenmc_%d", v)));
    vhrec_eff.emplace_back(hrec_eff);
    vhsel_eff.emplace_back(hsel_eff);
    vhrec_mat.emplace_back(hrec_mat);
    vhsel_mat.emplace_back(hsel_mat);
    vh2d.emplace_back(h2d);
  }
  Long64_t nentries = trec->GetEntries();
  for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
    trec->GetEntry(ientry);
    double e0_TeV = log10(e0 / 1000.);
    int i_type = TypeCor2TypeId(type);
    int ih2dbin = vh2d[i_type]->FindFixBin(e0_TeV, zenmc);
    double AtoN = DivideNoNan(kArea, vh2d[i_type]->GetBinContent(ih2dbin));
    vhrec_eff[i_type].Fill(e0_TeV, zenmc, AtoN);
    vhrec_mat[i_type].Fill(e0_TeV, zenmc, nfitc, AtoN);
    if (nfitc > 100 && compactness > 15) { // event selection
      vhsel_eff[i_type].Fill(e0_TeV, zenmc, AtoN);
      vhsel_mat[i_type].Fill(e0_TeV, zenmc, nfitc, AtoN);
    }
  }
  feff->Write();
  feff->Close();
}