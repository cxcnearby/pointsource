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

// RA = 83.633212;
// DEC = 22.014460;  //Crab

int main(int argc, char *argv[]) {

  if (argc < 2) {
    printf("%s  SourceType  [source RA]  [source DEC]\n(SourceType: 1 "
           "for Moon, 2 for Sun, 0 need RA & DEC)\n",
           argv[0]);
    exit(0);
  }
  int iSourceType = atoi(argv[1]);

  double dRA, dDEC, dZENITH, dAZI;

  if ((argc == 4) && (iSourceType == 0)) {
    dRA = atof(argv[2]);
    dDEC = atof(argv[3]);
  }

  clock_t ctStart, ctFinish;
  ctStart = clock();

  int n_zen = int(180. / zen_bin_width);
  int n_azi = int(360. / zen_bin_width);
  TFile *f = new TFile("crab_zen_dist.root", "recreate");
  TH1F *hzen = new TH1F("hzen", "zen dist", n_zen, 0, 180);
  TH1F *hazi = new TH1F("hazi", "azi dist", n_azi, -180, 180);

  iauASTROM astrom;
  double eo;
  double dut1 = -0.1641924;
  double xp = 152.182 * 0.001 * DAS2R;
  double yp = 275.050 * 0.001 * DAS2R;
  double mjd0 = 58800.;
  iauApco13(2400000.5, mjd0, dut1, 100.138794639 * D2R, 29.357656306 * D2R,
            4376.37, xp, yp, 0, 0, 0, 10, &astrom, &eo);
  for (int chana = 0; chana < 864000 * 365.2422 / 366.2422; ++chana) {
    // for (double chana = 0.; chana < 365.2422/366.2422; chana = chana + 1e-7)
    // {
    // use sidereal day instead of solar day to avoid overlap of target's track
    // on celestial sphere.
    double mjd = mjd0 + chana / 864000.;
    iauAper13(2400000.5, mjd, &astrom);
    switch (iSourceType) {
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
    // dAZI = dAZI * R2D;
    dZENITH = dZENITH * R2D;
    hzen->Fill(dZENITH, 0.1);
    hazi->Fill(dAZI, 0.1);
  }
  f->Write();
  f->Close();
  ctFinish = clock();
  cout << "use " << double((ctFinish - ctStart) / CLOCKS_PER_SEC) << " s."
       << endl;
  return 0;
}
