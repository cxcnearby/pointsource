/** @file evtselect.cc
 *  @author  changxc
 *  @date    2019.06.01
/* ========================================================= */
#include "TBranch.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TROOT.h"
#include "TTree.h"
#include "sofa.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

int main(int argc, char *argv[]) {

  if (argc < 5) {
    printf("%s  InputFileList  OutputFile.root  nfitc_min  compactness_min\n",
           argv[0]);
    exit(0);
  }
  string sInputFileList = argv[1];
  string sOutputFile = argv[2];
  const int kNFitcMin = atoi(argv[3]);
  const double kCompMin = atof(argv[4]);
  string::size_type pos = sOutputFile.rfind(".root");
  string sErrorFileList = sOutputFile.substr(0, pos) + ".error";
  ofstream ouError(sErrorFileList);

  clock_t ctStart, ctFinish;
  ctStart = clock();

  TString treename = "trec";
  TChain *cInput = new TChain(treename);
  string sXROOTD = "root://eos01.ihep.ac.cn/";
  bool emptylist = true;
  if (sInputFileList.find(".root") != string::npos) {
    cInput->Add(sInputFileList.c_str());
    emptylist = false;
  } else {
    string sFileName;
    string sServerFileName;
    ifstream inFile(sInputFileList);
    while (getline(inFile, sFileName)) {
      sServerFileName = sXROOTD + sFileName;
      TFile *testfile = TFile::Open(sServerFileName.c_str(), "READ");
      if (testfile && testfile->Get(treename)) {
        cInput->Add(sServerFileName.c_str());
        emptylist = false;
      } else {
        ouError << sFileName << endl;
      }
      if (testfile)
        testfile->Close();
    }
  }
  if (emptylist)
    exit(0);
  int nfitc;
  double compactness, recflag, chi2c, omega;
  cInput->SetBranchAddress("nfitc", &nfitc);
  cInput->SetBranchAddress("compactness", &compactness);
  cInput->SetBranchAddress("recflag", &recflag);
  cInput->SetBranchAddress("chi2c", &chi2c);
  cInput->SetBranchAddress("omega", &omega);
  TFile *fSelected = TFile::Open(sOutputFile.c_str(), "recreate");
  TTree *trec = cInput->CloneTree(0);
  Long64_t ntot = 0;
  Long64_t nsel = 0;
  TTree *tstat = new TTree("tstat", "some statistics");
  tstat->Branch("ntot", &ntot);
  tstat->Branch("nsel", &nsel);

  Long64_t nentries = cInput->GetEntries();
  ntot = nentries;
  int iFail;
  for (Long64_t i = 0; i < nentries; ++i) {
    iFail = cInput->GetEntry(i);
    if (iFail <= 0)
      continue;
    if (nfitc > kNFitcMin && compactness > kCompMin && recflag > 0 &&
        chi2c < 5 && omega < 5) {
      trec->Fill();
      ++nsel;
    }
  }
  delete cInput;
  tstat->Fill();
  trec->Write("", TObject::kOverwrite);
  tstat->Write("", TObject::kOverwrite);
  fSelected->Close();
  ctFinish = clock();
  cout << "total events: " << ntot << "; selected events: " << nsel
       << "; selecting use " << double((ctFinish - ctStart) / CLOCKS_PER_SEC)
       << " s." << endl;
  return 0;
}
