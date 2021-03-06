#include "constant.h"
#include "functions.h"

#include "TFile.h"
#include "TH2F.h"
#include "TTree.h"
#include <iostream>

using namespace std;

vector<vector<vector<double>>> SetStatMatrix(TString inputroot);
vector<vector<double>> SetRatio(const vector<double> &indexarray,
                                const int kNESmallbin);
void CalcEffareaResponse(TString inputroot,
                         const vector<vector<vector<double>>> &matrix,
                         const vector<vector<double>> &fluxratio);
inline double DivideNoNan(double up, double down) {
  return ((down < 1.e-8) ? 0 : up / down);
};

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << argv[0] << "  weight_added_root_file" << endl;
    exit(0);
  }
  TString input = argv[1];
  auto num_matrix = SetStatMatrix(input);
  auto suberange_ratio = SetRatio(kTypePowerIndex, kNESmallbin);
  CalcEffareaResponse(input, num_matrix, suberange_ratio);
  return 0;
}

vector<vector<vector<double>>> SetStatMatrix(TString inputroot) {
  vector<vector<vector<double>>> matrix(
      kTypeCode.size(),
      vector<vector<double>>(kEnergyBounds.size() - 1, vector<double>(2, 0.)));
  TFile *fsim = TFile::Open(inputroot, "read");
  TTree *tstat = (TTree *)fsim->Get("tstat");
  int type, erange;
  Long64_t sum_ntot, sum_nsel;
  tstat->SetBranchAddress("type", &type);
  tstat->SetBranchAddress("erange", &erange);
  tstat->SetBranchAddress("sum_ntot", &sum_ntot);
  tstat->SetBranchAddress("sum_nsel", &sum_nsel);
  Long64_t nentries = tstat->GetEntries();
  for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
    tstat->GetEntry(ientry);
    int i_erange = erange + 2;
    int i_type = TypeCor2TypeId(type);
    matrix[i_type][i_erange][0] = sum_ntot;
    matrix[i_type][i_erange][1] = sum_nsel;
  }
  fsim->Close();
  return matrix;
}

vector<vector<double>> SetRatio(const vector<double> &indexarray,
                                const int nsmallbin) {
  vector<vector<double>> ratio;
  const double step = pow(10., 1. / nsmallbin);
  for (unsigned int i = 0; i < indexarray.size(); ++i) {
    vector<double> type_ratio;
    double e1 = 1.;
    double inte = pow(e1, 1 + indexarray[i]) - pow(10, 1 + indexarray[i]);
    for (int j = 0; j < nsmallbin; ++j) {
      double bin_ratio =
          (pow(e1, 1 + indexarray[i]) - pow(e1 * step, 1 + indexarray[i])) /
          inte;
      type_ratio.emplace_back(bin_ratio);
      e1 *= step;
    }
    ratio.emplace_back(type_ratio);
  }
  return ratio;
}

void CalcEffareaResponse(TString inputroot,
                         const vector<vector<vector<double>>> &matrix,
                         const vector<vector<double>> &fluxratio) {
  TFile *fsim = TFile::Open(inputroot, "read");
  TTree *trec = (TTree *)fsim->Get("trec");
  float e0, zenmc, compactness, weight;
  int nfitc, type;
  trec->SetBranchAddress("e0", &e0);
  trec->SetBranchAddress("zenmc", &zenmc);
  trec->SetBranchAddress("compactness", &compactness);
  trec->SetBranchAddress("weight", &weight);
  trec->SetBranchAddress("nfitc", &nfitc);
  trec->SetBranchAddress("type", &type);
  TString outroot = inputroot.ReplaceAll(".root", "_effmat.root");
  TFile *feff = TFile::Open(outroot, "RECREATE");
  vector<TH2F> vhrec_eff, vhsel_eff, vhrec_mat, vhrec_mat_weight, vhsel_mat,
      vhsel_mat_weight;
  vector<TH1F> vhrec_num, vhrec_num_weight, vhsel_num, vhsel_num_weight;
  for (auto &&v : kTypeCode) {
    TH2F hrec_eff(Form("receff%d", v), Form("receff%d--energy x zenith", v),
                  kNEnergyBin, -2, 3, kNZenBin, 0,
                  kZenRange); // effarea of reconstructed events
    TH2F hsel_eff(Form("seleff%d", v), Form("seleff%d--energy x zenith", v),
                  kNEnergyBin, -2, 3, kNZenBin, 0,
                  kZenRange); // effarea of events after selection
    TH1F hrec_num(Form("recnum%d", v), Form("recnum%d--energy", v), kNEnergyBin,
                  -2, 3); // number of reconstructed events
    TH2F hrec_mat(Form("recmat%d", v), Form("recmat%d--energy x nfitc", v),
                  kNEnergyBin, -2, 3, kNPmtBin, 0,
                  kPmtRange); // response of reconstructed events
    TH1F hsel_num(Form("selnum%d", v), Form("selnum%d--energy", v), kNEnergyBin,
                  -2, 3); // number of events after selection
    TH2F hsel_mat(Form("selmat%d", v), Form("selmat%d--energy x nfitc", v),
                  kNEnergyBin, -2, 3, kNPmtBin, 0,
                  kPmtRange); // response of events after selection
    TH1F hrec_num_weight(Form("wrecnum%d", v), Form("wrecnum%d--energy", v),
                         kNEnergyBin, -2,
                         3); // weighted number of reconstructed events
    TH2F hrec_mat_weight(
        Form("wrecmat%d", v), Form("wrecmat%d--energy x nfitc", v), kNEnergyBin,
        -2, 3, kNPmtBin, 0,
        kPmtRange); // weighted response of reconstructed events
    TH1F hsel_num_weight(Form("wselnum%d", v), Form("wselnum%d--energy", v),
                         kNEnergyBin, -2,
                         3); // weighted number of events after selection
    TH2F hsel_mat_weight(
        Form("wselmat%d", v), Form("wselmat%d--energy x nfitc", v), kNEnergyBin,
        -2, 3, kNPmtBin, 0,
        kPmtRange); // weighted response of events after selection
    vhrec_eff.emplace_back(hrec_eff);
    vhsel_eff.emplace_back(hsel_eff);
    vhrec_num.emplace_back(hrec_num);
    vhrec_mat.emplace_back(hrec_mat);
    vhsel_num.emplace_back(hsel_num);
    vhsel_mat.emplace_back(hsel_mat);
    vhrec_num_weight.emplace_back(hrec_num_weight);
    vhrec_mat_weight.emplace_back(hrec_mat_weight);
    vhsel_num_weight.emplace_back(hsel_num_weight);
    vhsel_mat_weight.emplace_back(hsel_mat_weight);
  }
  vector<double> strip = StripArea(0., 180., kZenBinWidth);
  double totalarea = 2. * PI * (1. - cos(60. * D2R));
  Long64_t nentries = trec->GetEntries();
  for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
    trec->GetEntry(ientry);
    double e0_TeV = log10(e0 / 1000.);
    int i_erange = floor(e0_TeV) + 2;
    int i_sub_erange = floor((e0_TeV - floor(e0_TeV)) * kNESmallbin);
    int i_zen = int(zenmc / kZenBinWidth);
    int i_type = TypeCor2TypeId(type);
    double sum_ntot = matrix[i_type][i_erange][0];
    double sum_nsel = matrix[i_type][i_erange][1];
    double N_strip_tot =
        sum_ntot * fluxratio[i_type][i_sub_erange] * strip[i_zen] / totalarea;
    double AtoN = DivideNoNan(kArea, N_strip_tot);
    vhrec_eff[i_type].Fill(e0_TeV, zenmc, AtoN);
    vhrec_num[i_type].Fill(e0_TeV, 1);
    vhrec_mat[i_type].Fill(e0_TeV, nfitc, 1);
    vhrec_num_weight[i_type].Fill(e0_TeV, weight);
    vhrec_mat_weight[i_type].Fill(e0_TeV, nfitc, weight);
    if (nfitc > 100 && compactness > 15) { // event selection
      vhsel_eff[i_type].Fill(e0_TeV, zenmc, AtoN);
      vhsel_num[i_type].Fill(e0_TeV, 1);
      vhsel_mat[i_type].Fill(e0_TeV, nfitc, 1);
      vhsel_num_weight[i_type].Fill(e0_TeV, weight);
      vhsel_mat_weight[i_type].Fill(e0_TeV, nfitc, weight);
    }
  }
  fsim->Close();
  for (unsigned int i_type = 0; i_type < kTypeCode.size(); ++i_type) {
    for (int i = 1; i <= kNEnergyBin; ++i) {
      float rec_num = vhrec_num[i_type].GetBinContent(i);
      float sel_num = vhsel_num[i_type].GetBinContent(i);
      float rec_num_weight = vhrec_num_weight[i_type].GetBinContent(i);
      float sel_num_weight = vhsel_num_weight[i_type].GetBinContent(i);
      for (int j = 1; j <= kNPmtBin; ++j) {
        float rec_mat = vhrec_mat[i_type].GetBinContent(i, j);
        rec_mat = DivideNoNan(rec_mat, rec_num);
        vhrec_mat[i_type].SetBinContent(i, j, rec_mat);
        float sel_mat = vhsel_mat[i_type].GetBinContent(i, j);
        sel_mat = DivideNoNan(sel_mat, sel_num);
        vhsel_mat[i_type].SetBinContent(i, j, sel_mat);
        float rec_mat_weight = vhrec_mat_weight[i_type].GetBinContent(i, j);
        rec_mat_weight = DivideNoNan(rec_mat_weight, rec_num_weight);
        vhrec_mat_weight[i_type].SetBinContent(i, j, rec_mat_weight);
        float sel_mat_weight = vhsel_mat_weight[i_type].GetBinContent(i, j);
        sel_mat_weight = DivideNoNan(sel_mat_weight, sel_num_weight);
        vhsel_mat_weight[i_type].SetBinContent(i, j, sel_mat_weight);
      }
    }
  }
  feff->Write();
  feff->Close();
}
