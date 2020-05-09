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
#include "rebase/CopyTreeStruct.h"
#include "sofa.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

struct UsedLeaf {
  Int_t *I;
  Short_t *S;
  Long64_t *L;
  Float_t *F;
  Double_t *D;
  Bool_t *O;
}; // use pointer !!!

int main(int argc, char *argv[]) {

  if (argc < 5) {
    printf("%s  InputFileList  OutputFile.root  nfitc_min  "
           "compactness_min  [-t template.root]\n"
           "(-t template.root: -t need a template file name, "
           "or use the first-read root file as template)\n",
           argv[0]);
    exit(0);
  }
  string sInputFileList = argv[1];
  string sOutputFile = argv[2];
  const int kNFitcMin = atoi(argv[3]);
  const double kCompMin = atof(argv[4]);
  string template_name;
  for (int i = 5; i < argc; ++i) {
    if (strcmp(argv[i], "-t") == 0) {
      if (i + 2 <= argc) {
        template_name = argv[++i];
        TFile *testfile = TFile::Open(template_name.c_str(), "READ");
        if (!testfile) {
          printf("cannot open template file!\n");
          exit(0);
        }
        testfile->Close();
      } else {
        printf("need template.root name!\n");
        exit(0);
      }
    }
  }
  string::size_type pos = sOutputFile.rfind(".root");
  string sErrorFileList = sOutputFile.substr(0, pos) + ".error";
  ofstream ouError(sErrorFileList);

  Long64_t ntot = 0;
  Long64_t nsel = 0;
  clock_t ctStart, ctFinish;

  ctStart = clock();

  bool emptylist = true;
  CopyTreeStruct *cp;
  TChain *cInput;
  string sXROOTD = "root://eos01.ihep.ac.cn/";
  if (sInputFileList.find(".root") != string::npos) {
    if (template_name.empty()) {
      template_name = sInputFileList;
    }
    cp = new CopyTreeStruct(template_name.c_str(), "trec");
    int rSetReadTree = cp->SetReadTree(cInput);
    if (rSetReadTree) {
      std::cerr << "SetReadTree failed! error code: " << rSetReadTree
                << std::endl;
      exit(1);
    }
    cInput->Add(sInputFileList.c_str());
    emptylist = false;
  } else {
    string sFileName;
    string sServerFileName;
    ifstream inFile(sInputFileList);
    bool notemplate = true;
    while (getline(inFile, sFileName)) {
      sServerFileName = sXROOTD + sFileName;
      TFile *testfile = TFile::Open(sServerFileName.c_str(), "READ");
      if (template_name.empty() && testfile) {
        template_name = sServerFileName;
      }
      if (notemplate && !template_name.empty()) {
        cp = new CopyTreeStruct(template_name.c_str(), "trec");
        int rSetReadTree = cp->SetReadTree(cInput);
        if (rSetReadTree) {
          std::cerr << "SetReadTree failed! error code: " << rSetReadTree
                    << std::endl;
          exit(1);
        }
        notemplate = false;
      }
      if (testfile && testfile->Get(cp->GetTreeName())) {
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
  TFile *fSelected = TFile::Open(sOutputFile.c_str(), "recreate");
  TTree *tstat = new TTree("tstat", "some statistics");
  tstat->Branch("ntot", &ntot);
  tstat->Branch("nsel", &nsel);
  TTree *trec;
  int rSetWriteTree = cp->SetWriteTree(trec);
  if (rSetWriteTree) {
    std::cerr << "SetWriteTree failed! error code: " << rSetWriteTree
              << std::endl;
    exit(1);
  }
  UsedLeaf nfitc, compactness, recflag, chi2c, omega;
  cp->SetBond("nfitc", nfitc);
  cp->SetBond("compactness", compactness);
  cp->SetBond("recflag", recflag);
  cp->SetBond("chi2c", chi2c);
  cp->SetBond("omega", omega);
  Long64_t nentries = cInput->GetEntries();
  ntot = nentries;
  int iFail;
  for (Long64_t i = 0; i < nentries; ++i) {
    iFail = cInput->GetEntry(i);
    if (iFail <= 0)
      continue;
    if ((!nfitc.I || nfitc.I[0] > kNFitcMin) &&
        (!compactness.D || compactness.D[0] > kCompMin) &&
        (!recflag.D || recflag.D[0] > 0) && (!chi2c.D || chi2c.D[0] < 5) &&
        (!omega.D || omega.D[0] < 5)) {
      trec->Fill();
      ++nsel;
    }
  }
  delete cInput;
  tstat->Fill();
  tstat->Write("", TObject::kOverwrite);
  trec->Write("", TObject::kOverwrite);
  fSelected->Close();
  ctFinish = clock();
  cout << "total events: " << ntot << "; selected events: " << nsel
       << "; selecting use " << double((ctFinish - ctStart) / CLOCKS_PER_SEC)
       << " s." << endl;
  return 0;
}
