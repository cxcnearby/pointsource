/** @file evtselect.cc
 *  @author  changxc
 *  @date    2019.06.01
 * =========================================================
 */
#include "trackstat.h"
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
#include <functional>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
using namespace std;

void GenerateTrackStat(TString output,
                       const function<void(double, double *, double *)> &func);
function<void(double, double *, double *)> TargetEquaPosition(int sourcetype);

// RA = 83.633212;
// DEC = 22.014460;  //Crab
double dRA, dDEC, dZEN, dAZI;

int main(int argc, char *argv[]) {
  ArgInit arg(argc, argv);
  clock_t ctStart, ctFinish;
  ctStart = clock();
  auto func_pos = TargetEquaPosition(arg.source_type);
  GenerateTrackStat(arg.trackstatroot, func_pos);
  ctFinish = clock();
  cout << "use " << double((ctFinish - ctStart) / CLOCKS_PER_SEC) << " s."
       << endl;
  return 0;
}

void GenerateTrackStat(TString output,
                       const function<void(double, double *, double *)> &func) {
  int n_zen = int(180. / kZenBinWidth);
  int n_azi = int(360. / kZenBinWidth);
  TFile *f = new TFile(output, "recreate");
  TH1F *hzen = new TH1F("hzen", "zen dist", n_zen, 0, 180);
  TH1F *hazi = new TH1F("hazi", "azi dist", n_azi, -180, 180);
  iauASTROM astrom;
  double eo;
  double dut1 = -0.1641924;
  double xp = 152.182 * 0.001 * DAS2R;
  double yp = 275.050 * 0.001 * DAS2R;
  double mjd0 = 58800.; // 2019.11.13
  iauApco13(2400000.5, mjd0, dut1, 100.138794639 * D2R, 29.357656306 * D2R,
            4376.37, xp, yp, 0, 0, 0, 10, &astrom, &eo);
  for (int chana = 0; chana < 864000 * 365.2422 / 366.2422; ++chana) {
    // for (double chana = 0.; chana < 365.2422/366.2422; chana = chana + 1e-7)
    // {
    // use sidereal day instead of solar day to avoid overlap of target's track
    // on celestial sphere.
    double mjd = mjd0 + chana / 864000.;
    iauAper13(2400000.5, mjd, &astrom);
    func(mjd, &dRA, &dDEC);
    double ri, di, hob, dob, rob;
    iauAtciqz(dRA * D2R, dDEC * D2R, &astrom, &ri, &di);
    iauAtioq(ri, di, &astrom, &dAZI, &dZEN, &hob, &dob, &rob);
    dAZI = azirange(dAZI * R2D);
    // dAZI = dAZI * R2D;
    dZEN = dZEN * R2D;
    if (dZEN > 45.)
      continue;
    hzen->Fill(dZEN, 0.1);
    hazi->Fill(dAZI, 0.1);
  }
  f->Write();
  f->Close();
}

function<void(double, double *, double *)> TargetEquaPosition(int sourcetype) {
  switch (sourcetype) {
  case 1:
    return [](double mjd, double *ra, double *dec) {
      return moon_orbit(mjd, ra, dec);
    };
  case 2:
    return [](double mjd, double *ra, double *dec) {
      return sun_orbit(mjd, ra, dec);
    };
  default:
    return [](double mjd, double *ra, double *dec) {};
  }
}