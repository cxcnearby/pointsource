/** @file evtselect.cc
 *  @author  changxc
 *  @date    2019.06.01
 * =========================================================
 */
#include "Astro.h"
#include "constant.h"
#include "rebase/CopyTreeStruct.h"

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

struct UsedLeaf {
  Int_t *I;
  Short_t *S;
  Long64_t *L;
  Float_t *F;
  Double_t *D;
  Bool_t *O;
}; // use pointer !!!

struct WindowInfo {
  Int_t id;
  Double_t theta;
  Double_t phi;
  Double_t xd;
  Double_t yd;
};

int main(int argc, char *argv[]) {
  if (argc < 7) {
    printf("%s  InputFileList  OutputFile.root  WindowRadius  bgWindowNumber  "
           "ZenithMax  SourceType  [source RA]  [source DEC]  [-t "
           "template.root]\n"
           "(SourceType: -m for Moon, -s for Sun, -p need RA & DEC)\n"
           "(-t template.root: -t need a template file name, or use the "
           "first-read root file as template)\n",
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
  double dRA, dDEC, dZEN, dAZI;
  int source_type = 0;
  string template_name;
  for (int i = 6; i < argc; ++i) {
    if (strcmp(argv[i], "-m") == 0) {
      source_type = 1;
    } else if (strcmp(argv[i], "-s") == 0) {
      source_type = 2;
    } else if (strcmp(argv[i], "-p") == 0) {
      source_type = 3;
      if (i + 3 <= argc) {
        dRA = atof(argv[++i]);
        dDEC = atof(argv[++i]);
      } else {
        printf("wrong RA, DEC!\n");
        exit(0);
      }
    } else if (strcmp(argv[i], "-t") == 0) {
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
    } else {
      printf("wrong parameters!\n");
      exit(0);
    }
  }
  if (!source_type) {
    printf("source error!\n");
    exit(0);
  }
  FILE *fpLog = NULL;
  string::size_type pos = sOutputFile.rfind(".root");
  string sFpLog = sOutputFile.substr(0, pos) + ".log";
  if ((fpLog = fopen(sFpLog.c_str(), "w")) == NULL)
    printf("cannot open log file\n");
  string sErrorFileList = sOutputFile.substr(0, pos) + ".error";
  ofstream ouError(sErrorFileList);

  Long64_t ntot = 0;
  Long64_t nsel = 0;
  clock_t ctStart, ctFinish;
  double livetime = 0.;
  double expotime = 0.;

  ctStart = clock();

  bool emptylist = true;
  CopyTreeStruct *cp = nullptr;
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
  const int kNTimeBin = 86400 * 365;
  const double kDTime = 1. / 86400;
  TH1F *hlive = new TH1F("hlive", "events time dist. 10 second per bin",
                         kNTimeBin, 0, 365);
  TH1F *hexpo = new TH1F("hexpo", "events time dist. 10 second day per bin",
                         kNTimeBin, 0, 365);
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
  WindowInfo win;
  trec->Branch("id", &win.id);
  trec->Branch("theta", &win.theta);
  trec->Branch("phi", &win.phi);
  trec->Branch("xd", &win.xd);
  trec->Branch("yd", &win.yd);

  UsedLeaf mjd, zenc, azic;
  cp->SetBond("mjd", mjd);
  cp->SetBond("zenc", zenc);
  cp->SetBond("azic", azic);

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
      time_0 = int(mjd.D[0]);
    if (i == 0)
      iauApco13(2400000.5, mjd.D[0], dut1, 100.138794639 * D2R,
                29.357656306 * D2R, 4376.37, xp, yp, 0, 0, 0, 10, &astrom, &eo);
    // if (i%10000 == 0) iauApco13(2400000.5, mjd.D[0], dut1,
    // 100.138794639*D2R, 29.357656306*D2R, 4376.37, xp, yp, 0, 0, 0,
    // 10, &astrom, &eo);
    iauAper13(2400000.5, mjd.D[0], &astrom);
    double time_1 = mjd.D[0] - time_0;
    hlive->Fill(time_1);
    switch (source_type) {
    case 1:
      moon_orbit(mjd.D[0], &dRA, &dDEC);
      break;
    case 2:
      sun_orbit(mjd.D[0], &dRA, &dDEC);
      break;
    default:
      break;
    }
    double ri, di, hob, dob, rob;
    iauAtciqz(dRA * D2R, dDEC * D2R, &astrom, &ri, &di);
    iauAtioq(ri, di, &astrom, &dAZI, &dZEN, &hob, &dob, &rob);
    dAZI = azirange(dAZI * R2D);
    dZEN = dZEN * R2D;
    if (dZEN > kZenithMin && dZEN < kZenithMax) {
      hexpo->Fill(time_1);
      azic.D[0] = 60.54 - azic.D[0];
      azic.D[0] = azirange(azic.D[0]);
      for (win.id = -kWindowNumberOneSide; win.id <= kWindowNumberOneSide;
           win.id++) {
        double dAziDiff =
            2 * asin(sin(kWindowRadius * D2R) / sin(dZEN * D2R)) * R2D * win.id;
        double dAZIid = dAZI + dAziDiff;
        dAZIid = azirange(dAZIid);
        win.theta = distance_horizontal(dZEN, dAZIid, zenc.D[0], azic.D[0]);
        if (win.theta < kWindowRadius) {
          double daziid = azic.D[0] - dAziDiff;
          daziid = azirange(daziid);
          if (daziid < 0)
            daziid += 360.;
          double dra1, ddec1, dra, ddec;
          iauAtoiq("a", daziid * D2R, zenc.D[0] * D2R, &astrom, &dra1, &ddec1);
          iauAticq(dra1, ddec1, &astrom, &dra, &ddec);
          win.phi = direction_equatorial(dRA, dDEC, dra * R2D, ddec * R2D);
          win.xd = win.theta * sin(win.phi * D2R);
          win.yd = win.theta * cos(win.phi * D2R);
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
  tstat->Write("", TObject::kOverwrite);
  trec->Write("", TObject::kOverwrite);
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
