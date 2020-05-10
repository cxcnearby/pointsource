/** @file evtselect.cc
 *  @author  changxc
 *  @date    2019.06.01
/* ========================================================= */
#include "Astro.h"
#include "constant.h"

#include "sofa.h"

#include "TBranch.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TROOT.h"
#include "TTree.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

int main(int argc, char *argv[]) {
  if (argc < 7) {
    printf("%s  InputFileList  OutputFile.root  WindowRadius  bgWindowNumber  "
           "ZenithMax  SourceType  [source RA]  [source DEC]\n"
           "(SourceType: -m for Moon, -s for Sun, -p need RA & DEC)\n",
           argv[0]);
    exit(0);
  }
  string sInputFileList = argv[1];
  string sOutputFile = argv[2];
  const double kWindowRadius = atof(argv[3]);
  const int kWindowNumberOneSide = atof(argv[4]);
  const double kZenithMax = atof(argv[5]);
  const double kZenithMin = asin(sin(kWindowRadius * D2R) /
                                 sin(PI / (2 * kWindowNumberOneSide + 1))) *
                            R2D;
  double dRA, dDEC, dZENITH, dAZI;
  int source_type;
  bool ifsource = false;
  for (int i = 6; i < argc; ++i) {
    if (strcmp(argv[i], "-m") == 0) {
      source_type = 1;
      ifsource = true;
    } else if (strcmp(argv[i], "-s") == 0) {
      source_type = 2;
      ifsource = true;
    } else if (strcmp(argv[i], "-p") == 0) {
      source_type = 0;
      if (i + 3 <= argc) {
        dRA = atof(argv[++i]);
        dDEC = atof(argv[++i]);
        ifsource = true;
      } else {
        printf("wrong RA, DEC!\n");
        exit(0);
      }
    } else {
      printf("wrong parameters!\n");
      exit(0);
    }
  }
  // if (!ifsource) {
  //   printf("source error!\n");
  //   exit(0);
  // }
  FILE *fpLog = NULL;
  string::size_type pos = sOutputFile.rfind(".root");
  string sFpLog = sOutputFile.substr(0, pos) + ".log";
  if ((fpLog = fopen(sFpLog.c_str(), "w")) == NULL)
    printf("cannot open log file\n");
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
  double mjd, zenc, azic;
  cInput->SetBranchAddress("mjd", &mjd);
  cInput->SetBranchAddress("zenc", &zenc);
  cInput->SetBranchAddress("azic", &azic);
  TFile *fSelected = TFile::Open(sOutputFile.c_str(), "recreate");
  int id;
  double theta, phi, xd, yd;
  TTree *trec = cInput->CloneTree(0);
  trec->Branch("id", &id);
  trec->Branch("theta", &theta);
  trec->Branch("phi", &phi);
  trec->Branch("xd", &xd);
  trec->Branch("yd", &yd);
  Long64_t ntot = 0;
  Long64_t nsel = 0;
  TTree *tstat = new TTree("tstat", "some statistics");
  tstat->Branch("ntot", &ntot);
  tstat->Branch("nsel", &nsel);
  const int kNTimeBin = 86400 * 365;
  const double kDTime = 1. / 86400;
  double livetime = 0.;
  double expotime = 0.;
  TH1F *hlive = new TH1F("hlive", "events time dist. 10 second per bin",
                         kNTimeBin, 0, 365);
  TH1F *hexpo = new TH1F("hexpo", "events time dist. 10 second day per bin",
                         kNTimeBin, 0, 365);

  iauASTROM astrom;
  double eo;
  const double dut1 = -0.1641924;
  const double xp = 152.182 * 0.001 * DAS2R;
  const double yp = 275.050 * 0.001 * DAS2R;
  Long64_t nentries = cInput->GetEntries();
  ntot = nentries;
  int time_0;
  int iFail;
  for (Long64_t i = 0; i < nentries; ++i) {
    iFail = cInput->GetEntry(i);
    if (iFail <= 0)
      continue;
    if (i == 0)
      time_0 = int(mjd);
    if (i == 0)
      iauApco13(2400000.5, mjd, dut1, 100.138794639 * D2R, 29.357656306 * D2R,
                4376.37, xp, yp, 0, 0, 0, 10, &astrom, &eo);
    // if (i%10000 == 0) iauApco13(2400000.5, mjd, dut1,
    // 100.138794639*D2R, 29.357656306*D2R, 4376.37, xp, yp, 0, 0, 0,
    // 10, &astrom, &eo);
    iauAper13(2400000.5, mjd, &astrom);
    double time_1 = mjd - time_0;
    hlive->Fill(time_1);
    switch (source_type) {
    case 1:
      moon_orbit(mjd, &dRA, &dDEC);
      break;
    case 2:
      sun_orbit(mjd, &dRA, &dDEC);
      break;
    default:
      break;
    }
    double ri, di, hob, dob, rob;
    iauAtciqz(dRA * D2R, dDEC * D2R, &astrom, &ri, &di);
    iauAtioq(ri, di, &astrom, &dAZI, &dZENITH, &hob, &dob, &rob);
    dAZI = azirange(dAZI * R2D);
    dZENITH = dZENITH * R2D;
    if (dZENITH > kZenithMin && dZENITH < kZenithMax) {
      hexpo->Fill(time_1);
      azic = 60.54 - azic;
      azic = azirange(azic);
      for (id = -kWindowNumberOneSide; id <= kWindowNumberOneSide; id++) {
        double dAziDiff =
            2 * asin(sin(kWindowRadius * D2R) / sin(dZENITH * D2R)) * R2D * id;
        double dAZIid = dAZI + dAziDiff;
        dAZIid = azirange(dAZIid);
        theta = distance_horizontal(dZENITH, dAZIid, zenc, azic);
        if (theta < kWindowRadius) {
          double daziid = azic - dAziDiff;
          daziid = azirange(daziid);
          if (daziid < 0)
            daziid += 360.;
          double dra1, ddec1, dra, ddec;
          iauAtoiq("a", daziid * D2R, zenc * D2R, &astrom, &dra1, &ddec1);
          iauAticq(dra1, ddec1, &astrom, &dra, &ddec);
          phi = direction_equatorial(dRA, dDEC, dra * R2D, ddec * R2D);
          xd = theta * sin(phi * D2R);
          yd = theta * cos(phi * D2R);
          trec->Fill();
          ++nsel;
          break;
        }
      }
    }
  }
  for (int i = 1; i <= kNTimeBin; i++) {
    double n_live = hlive->GetBinContent(i);
    double n_expo = hexpo->GetBinContent(i);
    if (n_live > 0.1)
      livetime += kDTime;
    if (n_expo > 0.1)
      expotime += kDTime;
  }
  delete hlive;
  delete hexpo;
  delete cInput;
  tstat->Fill();
  trec->Write("", TObject::kOverwrite);
  tstat->Write("", TObject::kOverwrite);
  fSelected->Close();
  ctFinish = clock();
  fprintf(fpLog, "%f %f\n", livetime, expotime);
  if (fpLog)
    fclose(fpLog);
  cout << "LiveTime: " << livetime << " d; ExpoTime: " << expotime
       << " d; total events: " << ntot << "; selected events: " << nsel
       << "; selecting use " << double((ctFinish - ctStart) / CLOCKS_PER_SEC)
       << " s." << endl;
  return 0;
}
