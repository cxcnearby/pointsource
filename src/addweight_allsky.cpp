#include "constant.h"
#include "functions.h"

#include "TBranch.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2.h"
#include "TNtupleD.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TTree.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << argv[0] << "  input.root" << endl;
    exit(1);
  }
  string input = argv[1];
  string outroot = input.substr(0, input.length() - 5) + "wt.root";

  clock_t ctStart, ctFinish;
  ctStart = clock();

  TString treename = "trec";
  TFile *fin = TFile::Open(input.c_str(), "READ");
  TTree *t1 = (TTree *)fin->Get(treename);
  Int_t type, nfitc;
  Float_t zenmc, e0;
  t1->SetBranchAddress("type", &type);
  t1->SetBranchAddress("nfitc", &nfitc);
  t1->SetBranchAddress("zenmc", &zenmc);
  t1->SetBranchAddress("e0", &e0);
  t1->GetEntry(0);
  treename = "tstat";
  TTree *t2 = (TTree *)fin->Get(treename);
  Long64_t ntot, nsel;
  t2->SetBranchAddress("ntot", &ntot);
  t2->SetBranchAddress("nsel", &nsel);

  TFile *foroot = new TFile(outroot.c_str(), "recreate");
  TH1F *h_e0 = (TH1F *)fin->Get(Form("h_e0_%d", type))->Clone();
  TH2F *h_e0_zenmc = (TH2F *)fin->Get(Form("h_e0_zenmc_%d", type))->Clone();
  TH3F *h_e0_zenmc_nfitc =
      (TH3F *)fin->Get(Form("h_e0_zenmc_nfitc_%d", type))->Clone();
  foroot->WriteTObject(h_e0, "WriteDelete");
  foroot->WriteTObject(h_e0_zenmc, "WriteDelete");
  foroot->WriteTObject(h_e0_zenmc_nfitc, "WriteDelete");
  delete h_e0;
  delete h_e0_zenmc;
  delete h_e0_zenmc_nfitc;
  float weight;
  TTree *trec = t1->CloneTree(0);
  trec->Branch("weight", &weight);
  int erange;
  Long64_t sum_ntot = 0;
  Long64_t sum_nsel = 0;
  TTree *tstat = new TTree("tstat", "some statistics");
  tstat->Branch("type", &type);
  tstat->Branch("erange", &erange);
  tstat->Branch("sum_ntot", &sum_ntot);
  tstat->Branch("sum_nsel", &sum_nsel);
  Long64_t nentries;
  nentries = t2->GetEntries();
  for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
    t2->GetEntry(ientry);
    sum_ntot += ntot;
    sum_nsel += nsel;
  }
  delete t2;
  vector<double> strip = StripArea(0., 180., kZenBinWidth);
  vector<double> bin_flux =
      BinnedIntegratedFlux(IntegratedFlux(Flux(type)), kEnergyBounds);
  double totalarea = 2. * PI * (1. - cos(60. * D2R));
  nentries = t1->GetEntries();
  for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
    t1->GetEntry(ientry);
    erange = floor(log10(e0 / 1000.));
    int i_e0 = erange + 2;
    int i_zen = int(zenmc / kZenBinWidth);
    if (type == 0) {
      weight = 0.;
    } else {
      weight = bin_flux[i_e0] * strip[i_zen] * 86400. * kArea *
               cos((i_zen + 0.5) * kZenBinWidth * D2R) /
               (sum_ntot * strip[i_zen] / totalarea);
    }
    trec->Fill();
  }
  delete t1;
  tstat->Fill();
  foroot->WriteTObject(trec, "WriteDelete");
  foroot->WriteTObject(tstat, "WriteDelete");
  foroot->Close();
  ctFinish = clock();
  printf("Adding branches use %d s\n",
         int((ctFinish - ctStart) / CLOCKS_PER_SEC));
  return 0;
}
