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

double corsika_event_number_rb[] = {1.e-6, 1.e8, 1.e8, 2.e6, 4.e5};
double corsika_event_number_rc[] = {1.e-6, 9.914e7, 9.9963e7, 1.9998e6,
                                    3.9988e5};

int main(int argc, char *argv[]) {
  if (argc < 5) {
    cout << argv[0]
         << "  input.root  particle_type  window_radius  direction_error"
         << endl;
    exit(1);
  }
  string input = argv[1];
  int type = atoi(argv[2]);
  const double window_radius = atof(argv[3]);
  const double direction_error = atof(argv[4]);
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
  TFile *fcrab = new TFile("crab_zen_dist.root", "read");
  TH1F *hzen = (TH1F *)fcrab->Get("hzen");

  TFile *foroot = new TFile(outroot.c_str(), "recreate");
  double weight;
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
  for (int ientry = 0; ientry < nentries; ++ientry) {
    cInput2->GetEntry(ientry);
    sum_ntot += ntot;
    sum_nsel += nsel;
  }
  delete cInput2;
  vector<double> strip = striparea();
  vector<double> bin_flux = binned_integrated_flux(type);
  vector<double> point_time_zen = point_duration_of_zenith_bin();
  vector<double> inwindow_time_zen =
      inwindow_duration_of_zenith_bin(window_radius, direction_error);
  double totalarea = 2. * PI * (1. - cos(60. * D2R));
  nentries = cInput1->GetEntries();
  for (int ientry = 0; ientry < nentries; ++ientry) {
    cInput1->GetEntry(ientry);
    erange = floor(log10(e0 / 1000.));
    int i_e0 = erange + 2;
    int i_zen = int(zenmc / kZenBinWidth);
    if (type == 0) {
      weight = bin_flux[i_e0] * point_time_zen[i_zen] * kArea *
               cos((i_zen + 0.5) * kZenBinWidth * D2R) /
               (sum_ntot * strip[i_zen] / totalarea);
    } else {
      weight = bin_flux[i_e0] * strip[i_zen] * inwindow_time_zen[i_zen] *
               kArea * cos((i_zen + 0.5) * kZenBinWidth * D2R) /
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
