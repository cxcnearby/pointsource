/** @file evtselect.cc
 *  @author  changxc
 *  @date    2019.06.01
 * =========================================================
 */
#include "Astro.h"
#include "constant.h"
#include "functions.h"

#include "sofa.h"

#include "TBranch.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TLeaf.h"
#include "TROOT.h"
#include "TTree.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

Int_t u[14];
Float_t v[302];

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
  ifstream inFile(sInputFileList);
  ofstream ouError(sErrorFileList);

  clock_t ctStart, ctFinish;
  ctStart = clock();

  TString treename = "t_evth";
  TChain *cInput = new TChain(treename);
  string sXROOTD = "root://eos01.ihep.ac.cn/";
  bool emptylist = true;
  string sFileName;
  string sServerFileName;
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
  if (emptylist)
    exit(0);
  cInput->SetBranchAddress("u", u);
  cInput->SetBranchAddress("v", v);
  TLeaf *leaf = cInput->GetLeaf("v");
  int v_size = leaf->GetNdata();
  const int kVersion = v_size < 300 ? 1 : 2;

  int type = GetTypeFromDir(sInputFileList);
  TFile *fSelected = TFile::Open(sOutputFile.c_str(), "recreate");
  TH1F *h_e0 = new TH1F(Form("h_e0_%d", type), Form("log10(e0) %d", type),
                        kNEnergyBin, -2, 3);
  TH2F *h_e0_zenmc =
      new TH2F(Form("h_e0_zenmc_%d", type), Form("log10(e0) x zenmc %d", type),
               kNEnergyBin, -2, 3, kNZenBin, 0, kZenRange);
  TH3F *h_e0_zenmc_nfitc =
      new TH3F(Form("h_e0_zenmc_nfitc_%d", type),
               Form("log10(e0) x zenmc x nfitc %d", type), kNEnergyBin, -2, 3,
               kNZenBin, 0, kZenRange, kNPmtBin, 0, kPmtRange);
  TTree *trec = new TTree("trec", "par selected events");
  trec->Branch("type", &type);
  trec->Branch("npmtinshower", &u[1]);
  trec->Branch("noise", &u[2]);
  trec->Branch("nhit", &u[3]);
  trec->Branch("nfitp", &u[4]);
  trec->Branch("nfitc", &u[5]);
  trec->Branch("ndetc", &u[6]);
  trec->Branch("npea", &u[7]);
  trec->Branch("npec", &u[8]);
  trec->Branch("multi", &u[9]);
  trec->Branch("recflag1", &u[11]);
  trec->Branch("recflag2", &u[12]);
  trec->Branch("e0", &v[3]);
  trec->Branch("zenmc", &v[10]);
  trec->Branch("azimc", &v[11]);
  trec->Branch("zen", &v[275]);
  trec->Branch("azi", &v[276]);
  trec->Branch("chi2p", &v[277]);
  trec->Branch("xc1", &v[278]);
  trec->Branch("yc1", &v[279]);
  trec->Branch("rmd20", &v[280]);
  trec->Branch("dcore", &v[281]);
  trec->Branch("zenc", &v[282]);
  trec->Branch("azic", &v[283]);
  trec->Branch("omega", &v[284]);
  trec->Branch("chi2c", &v[285]);
  Long64_t ntot = 0;
  Long64_t nsel = 0;
  TTree *tstat = new TTree("tstat", "some statistics");
  tstat->Branch("ntot", &ntot);
  tstat->Branch("nsel", &nsel);

  switch (kVersion) {
  case 1:
    trec->Branch("compactness", &v[286]);
    trec->Branch("pincness", &v[287]);
    // trec->Branch("compactness1", &v[288]);
    // trec->Branch("compactness2", &v[289]);
    // trec->Branch("compactness3", &v[290]);
    // trec->Branch("compactness4", &v[291]);
    break;
  case 2:
    trec->Branch("xc", &v[286]);
    trec->Branch("yc", &v[287]);
    trec->Branch("compactness", &v[288]);
    trec->Branch("compactness2", &v[289]);
    trec->Branch("compactness3", &v[290]);
    trec->Branch("compactness4", &v[291]);
    trec->Branch("pincness", &v[292]);
    break;
  default:
    break;
  }

  Long64_t nentries = cInput->GetEntries();
  ntot = nentries;
  int iFail;
  switch (kVersion) {
  case 1:
    for (Long64_t i = 0; i < nentries; ++i) {
      iFail = cInput->GetEntry(i);
      if (iFail <= 0)
        continue;
      if (u[5] > kNFitcMin && v[286] > kCompMin && u[11] > 0 && v[285] < 5) {
        v[10] *= R2D;
        v[11] = v[11] * R2D - 90.;
        // v[276] -= - 90.;  // TODO: these two azi seem have been converted.
        // v[283] -= 90.;
        v[278] += 75.5;
        v[279] += 56.5;
        h_e0->Fill(v[3]);
        h_e0_zenmc->Fill(v[3], v[10]);
        h_e0_zenmc_nfitc->Fill(v[3], v[10], u[5]);
        trec->Fill();
        ++nsel;
      }
    }
    break;
  case 2:
    for (Long64_t i = 0; i < nentries; ++i) {
      iFail = cInput->GetEntry(i);
      if (iFail <= 0)
        continue;
      if (u[5] > kNFitcMin && v[288] > kCompMin && u[11] > 0 && v[285] < 5) {
        v[10] *= R2D;
        v[11] = v[11] * R2D - 90.;
        v[278] += 75.5;
        v[279] += 56.5;
        v[286] += 75.5;
        v[287] += 56.5;
        h_e0->Fill(v[3]);
        h_e0_zenmc->Fill(v[3], v[10]);
        h_e0_zenmc_nfitc->Fill(v[3], v[10], u[5]);
        trec->Fill();
        ++nsel;
      }
    }
    break;
  default:
    break;
  }
  delete cInput;
  tstat->Fill();
  fSelected->WriteTObject(h_e0, "WriteDelete");
  fSelected->WriteTObject(h_e0_zenmc, "WriteDelete");
  fSelected->WriteTObject(h_e0_zenmc_nfitc, "WriteDelete");
  trec->Write("", TObject::kOverwrite);
  tstat->Write("", TObject::kOverwrite);
  fSelected->Close();
  ctFinish = clock();
  cout << "total events: " << ntot << "; selected events: " << nsel
       << "; selecting use " << double((ctFinish - ctStart) / CLOCKS_PER_SEC)
       << " s." << endl;
  return 0;
}
