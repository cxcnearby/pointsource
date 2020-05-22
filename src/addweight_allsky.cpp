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
  if (argc < 3) {
    cout << argv[0] << "  input.root  particle_type" << endl;
    exit(1);
  }
  string input = argv[1];
  int type = atoi(argv[2]);
  string outroot = input.substr(0, input.length() - 5) + "wt.root";

  clock_t ctStart, ctFinish;
  ctStart = clock();

  TString treename = "trec";
  TChain *cInput1 = new TChain(treename);
  TFile *testfile = TFile::Open(input.c_str(), "READ");
  if (testfile && testfile->Get(treename)) {
    testfile->Close();
    cInput1->Add(input.c_str());
  } else {
    cout << "wrong file!" << endl;
    exit(1);
  }
  Float_t zenmc, e0;
  cInput1->SetBranchAddress("zenmc", &zenmc);
  cInput1->SetBranchAddress("e0", &e0);
  treename = "tstat";
  TChain *cInput2 = new TChain(treename);
  cInput2->Add(input.c_str());
  Long64_t ntot, nsel;
  cInput2->SetBranchAddress("ntot", &ntot);
  cInput2->SetBranchAddress("nsel", &nsel);

  TFile *foroot = new TFile(outroot.c_str(), "recreate");
  float weight;
  TTree *trec = cInput1->CloneTree(0);
  trec->Branch("type", &type);
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
  nentries = cInput2->GetEntries();
  for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
    cInput2->GetEntry(ientry);
    sum_ntot += ntot;
    sum_nsel += nsel;
  }
  delete cInput2;
  vector<double> strip = StripArea(0., 180., kZenBinWidth);
  vector<double> bin_flux =
      BinnedIntegratedFlux(IntegratedFlux(Flux(type)), kEnergyBounds);
  double totalarea = 2. * PI * (1. - cos(60. * D2R));
  nentries = cInput1->GetEntries();
  for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
    cInput1->GetEntry(ientry);
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
  delete cInput1;
  tstat->Fill();
  foroot->cd();
  tstat->Write("", TObject::kOverwrite);
  trec->Write("", TObject::kOverwrite);
  foroot->Close();
  ctFinish = clock();
  printf("Adding branches use %d s\n",
         int((ctFinish - ctStart) / CLOCKS_PER_SEC));
  return 0;
}
