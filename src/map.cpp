#include "Astro.h"
#include "TBranch.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "constant.h"
#include <TH1.h>
#include <TH2.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include "ener_reso.h"

using namespace std;

double sig_smooth(vector<vector<float>> &vOn, vector<vector<float>> &vOff,
                  int x, int y, int id, double binwidth, int nbin,
                  double r_smooth, double *Non, double *Noff);
double gaus32(int kk, double r);
double bigaussian_mean(double theta, double sigma);
double bigaussian_simple(double theta, double sigma);

struct s_wcdaeventsSel {
  double mjd;
  int nhit;
  double npea;
  double zen;
  double azi;
  double xc;
  double yc;
  int ndetc;
  int nfitc;
  double npec;
  double zenc;
  double azic;
  double omega;
  double chi2p;
  double chi2c;
  double compactness;
  double pincness;
  int id;
  double theta;
  double phi;
  double xd;
  double yd;
};

int nk = 1;

int main(int argc, char *argv[]) {

  if (argc < 7) {
    printf("%s  InputFile.root  OutputFile.root  WindowRadius  bgWindowNumber  "
           "BinWidth  SmoothRadius\nUse \\* instead of * if appearing in "
           "InputFile.root\n",
           argv[0]);
    exit(0);
  }
  string sInputFile = argv[1];
  string sOutputFile = argv[2];
  double dWindowRadius = atof(argv[3]);
  int iWindowNumberOneSide = atof(argv[4]);
  double dBinWidth = atof(argv[5]);
  double dBinWidth1 = 0.1;
  double dSmoothRadius = atof(argv[6]);
  int iBinNumber = int(dWindowRadius / dBinWidth * 2.0 + 1.0);

  vector<vector<float>> vOn(iBinNumber + 1, vector<float>(iBinNumber + 1));
  vector<vector<float>> vOff(iBinNumber + 1, vector<float>(iBinNumber + 1));

  s_wcdaeventsSel s_EventSel;
  double dMapRange = 3.0;
  int iMapBinNumber = int(dMapRange / dBinWidth * 2.0 + 1.0);
  int iBinShift = (iBinNumber - iMapBinNumber) / 2;
  double dTotalOn = 0.0;
  double dTotalOff = 0.0;
  double dSig, dSigSmooth;
  double dNOn, dNOff;

  TChain *cInput = new TChain("trec");
  cInput->Add(sInputFile.c_str());

  cInput->SetBranchAddress("mjd", &s_EventSel.mjd);
  cInput->SetBranchAddress("nhit", &s_EventSel.nhit);
  cInput->SetBranchAddress("npea", &s_EventSel.npea);
  cInput->SetBranchAddress("zen", &s_EventSel.zen);
  cInput->SetBranchAddress("azi", &s_EventSel.azi);
  cInput->SetBranchAddress("xc1", &s_EventSel.xc);
  cInput->SetBranchAddress("yc1", &s_EventSel.yc);
  cInput->SetBranchAddress("ndetc", &s_EventSel.ndetc);
  cInput->SetBranchAddress("nfitc", &s_EventSel.nfitc);
  cInput->SetBranchAddress("npec", &s_EventSel.npec);
  cInput->SetBranchAddress("zenc", &s_EventSel.zenc);
  cInput->SetBranchAddress("azic", &s_EventSel.azic);
  cInput->SetBranchAddress("omega", &s_EventSel.omega);
  cInput->SetBranchAddress("chi2p", &s_EventSel.chi2p);
  cInput->SetBranchAddress("chi2c", &s_EventSel.chi2c);
  cInput->SetBranchAddress("compactness", &s_EventSel.compactness);
  cInput->SetBranchAddress("pincness", &s_EventSel.pincness);
  cInput->SetBranchAddress("id", &s_EventSel.id);
  cInput->SetBranchAddress("theta", &s_EventSel.theta);
  cInput->SetBranchAddress("phi", &s_EventSel.phi);
  cInput->SetBranchAddress("xd", &s_EventSel.xd);
  cInput->SetBranchAddress("yd", &s_EventSel.yd);

  TFile *fMap = new TFile(sOutputFile.c_str(), "recreate");
  TH2F *hOn = new TH2F(
      "hon", "on events", iBinNumber, -(dWindowRadius + dBinWidth / 2.0),
      dWindowRadius + dBinWidth / 2.0, iBinNumber,
      -(dWindowRadius + dBinWidth / 2.0), dWindowRadius + dBinWidth / 2.0);
  TH2F *hOff = new TH2F(
      "hoff", "off events", iBinNumber, -(dWindowRadius + dBinWidth / 2.0),
      dWindowRadius + dBinWidth / 2.0, iBinNumber,
      -(dWindowRadius + dBinWidth / 2.0), dWindowRadius + dBinWidth / 2.0);
  TH2F *hExc = new TH2F(
      "hexc", "on-off events", iBinNumber, -(dWindowRadius + dBinWidth / 2.0),
      dWindowRadius + dBinWidth / 2.0, iBinNumber,
      -(dWindowRadius + dBinWidth / 2.0), dWindowRadius + dBinWidth / 2.0);
  // TH2F *hExc = new TH2F(
  //     "hexc", "on-off events", iMapBinNumber, -(dMapRange + dBinWidth / 2.0),
  //     dMapRange + dBinWidth / 2.0, iMapBinNumber,
  //     -(dMapRange + dBinWidth / 2.0), dMapRange + dBinWidth / 2.0);
  TH1F *hSig1 = new TH1F("hsig1", "1-D significance", 300, -15.0, 15.0);
  TH2F *hOnS =
      new TH2F("hons", "on events (smoothed)", iMapBinNumber,
               -(dMapRange + dBinWidth / 2.0), dMapRange + dBinWidth / 2.0,
               iMapBinNumber, -(dMapRange + dBinWidth / 2.0),
               dMapRange + dBinWidth / 2.0);
  TH2F *hOffS =
      new TH2F("hoffs", "off events (smoothed)", iMapBinNumber,
               -(dMapRange + dBinWidth / 2.0), dMapRange + dBinWidth / 2.0,
               iMapBinNumber, -(dMapRange + dBinWidth / 2.0),
               dMapRange + dBinWidth / 2.0);
  TH2F *hExcS =
      new TH2F("hexcs", "on-off events (smoothed)", iMapBinNumber,
               -(dMapRange + dBinWidth / 2.0), dMapRange + dBinWidth / 2.0,
               iMapBinNumber, -(dMapRange + dBinWidth / 2.0),
               dMapRange + dBinWidth / 2.0);
  TH1F *hSigS1 =
      new TH1F("hsigs1", "1-D significance (smoothed)", 300, -15.0, 15.0);
  TH2F *hSigS2 =
      new TH2F("hsigs2", "2-D significance (smoothed)", iMapBinNumber,
               -(dMapRange + dBinWidth / 2.0), dMapRange + dBinWidth / 2.0,
               iMapBinNumber, -(dMapRange + dBinWidth / 2.0),
               dMapRange + dBinWidth / 2.0);
  TH1F *hRa = new TH1F("hra", "on-off events", 2 * dMapRange / dBinWidth1 + 1,
                       -(3 + dBinWidth1 / 2.0), 3 + dBinWidth1 / 2.0);
  TH1F *hDec = new TH1F("hdec", "on-off events", 2 * dMapRange / dBinWidth1 + 1,
                        -(3 + dBinWidth1 / 2.0), 3 + dBinWidth1 / 2.0);

  Long64_t nentries = cInput->GetEntries();
  for (Long64_t i = 0; i < nentries; i++) {
    cInput->GetEntry(i);
    if (s_EventSel.nfitc > 100 && s_EventSel.compactness > 15. &&
        sqrt(s_EventSel.xc * s_EventSel.xc + s_EventSel.yc * s_EventSel.yc) <
            400000) {
      if (s_EventSel.id == 0) {
        hOn->Fill(s_EventSel.xd, s_EventSel.yd);
        dTotalOn++;
      } else {
        hOff->Fill(s_EventSel.xd, s_EventSel.yd);
        dTotalOff++;
      }
    }
  }
  delete cInput;

  for (int i = 1; i <= iBinNumber; i++) {
    for (int j = 1; j <= iBinNumber; j++) {
      vOn[i][j] = hOn->GetBinContent(i, j);
      vOff[i][j] = hOff->GetBinContent(i, j);
      hExc->SetBinContent(i, j,
                          vOn[i][j] - vOff[i][j] / (2 * iWindowNumberOneSide));
    }
  }

  for (int i = 1; i <= iMapBinNumber; i++) {
    for (int j = 1; j <= iMapBinNumber; j++) {
      dSigSmooth = sig_smooth(vOn, vOff, i + iBinShift, j + iBinShift,
                              iWindowNumberOneSide, dBinWidth, iBinNumber,
                              dSmoothRadius, &dNOn, &dNOff);
      // hExc->SetBinContent(i, j,
      //                     vOn[i + iBinShift][j + iBinShift] - vOff[i +
      // iBinShift][j + iBinShift] / (2 * iWindowNumberOneSide));
      if (vOn[i + iBinShift][j + iBinShift] != 0)
        dSig =
            (vOn[i + iBinShift][j + iBinShift] -
             vOff[i + iBinShift][j + iBinShift] / (2 * iWindowNumberOneSide)) /
            sqrt(vOn[i + iBinShift][j + iBinShift] +
                 vOff[i + iBinShift][j + iBinShift] /
                     (2 * iWindowNumberOneSide * 2 * iWindowNumberOneSide));
      else
        dSig = -100;
      hSig1->Fill(dSig);
      hOnS->SetBinContent(i, j, dNOn);
      hOffS->SetBinContent(i, j, dNOff);
      hExcS->SetBinContent(i, j, dNOn - dNOff / (2 * iWindowNumberOneSide));
      hSigS1->Fill(dSigSmooth);
      hSigS2->SetBinContent(i, j, dSigSmooth);
      if (j >= (iMapBinNumber + 1) / 2 - 30 &&
          j <= (iMapBinNumber + 1) / 2 + 30)
        hRa->Fill((i - 1) * dBinWidth - dMapRange,
                  vOn[i + iBinShift][j + iBinShift] -
                      vOff[i + iBinShift][j + iBinShift] /
                          (2 * iWindowNumberOneSide));
      if (i >= (iMapBinNumber + 1) / 2 - 30 &&
          i <= (iMapBinNumber + 1) / 2 + 30)
        hDec->Fill((j - 1) * dBinWidth - dMapRange,
                   vOn[i + iBinShift][j + iBinShift] -
                       vOff[i + iBinShift][j + iBinShift] /
                           (2 * iWindowNumberOneSide));
      if (i == (iMapBinNumber + 1) / 2 && j == (iMapBinNumber + 1) / 2)
        printf("center %f %f %f %f %f %f %f %f %f\n", dSmoothRadius, dSigSmooth,
               dNOn, dNOff, dTotalOn, dTotalOff,
               (dNOn - dNOff / (2 * iWindowNumberOneSide)),
               2. * iWindowNumberOneSide,
               (dNOn - dNOff / (2 * iWindowNumberOneSide)) /
                   sqrt(dNOff / (2 * iWindowNumberOneSide)));
    }
  }
  fMap->Write();
  fMap->Close();
  return 0;
}

double sig_smooth(vector<vector<float>> &vOn, vector<vector<float>> &vOff,
                  int x, int y, int id, double binwidth, int nbin,
                  double r_smooth, double *Non, double *Noff) {
  int n_r;
  int i, j, i0, j0, i1, j1;
  double w = 0, r = 0;
  double dN = 0, dNoff = 0;
  double dis_square;
  double r_smooth_square = r_smooth * r_smooth;
  n_r = int(1. * r_smooth / binwidth);

  *Non = 0.;
  *Noff = 0.;

  i0 = x - n_r;
  i1 = x + n_r;
  j0 = y - n_r;
  j1 = y + n_r;

  if (i0 < 1)
    i0 = 1;
  if (i1 > nbin)
    i1 = nbin;

  if (j0 < 1)
    j0 = 1;
  if (j1 > nbin)
    j1 = nbin;

  for (i = i0; i <= i1; i++) {
    for (j = j0; j <= j1; j++) {

      dis_square =
          ((i - x) * (i - x) + (j - y) * (j - y)) * binwidth * binwidth;
      if (dis_square > 1.5 * r_smooth_square)
        continue;

      if (nk > 12)
        w = 1;
      else
        // w = gaus32(nk, sqrt(dis_square) * D2R);
        w = 1 * pow(bigaussian_simple(sqrt(dis_square), r_smooth), 2);
      *Non += vOn[i][j] * w;
      *Noff += vOff[i][j] * w;
      dN += vOn[i][j] * w * w;
      dNoff += vOff[i][j] * w * w;
    }
  }

  if (dN != 0)
    // return (*Non - *Noff / (2 * id)) / sqrt(dN + dNoff / (2 * 2 * id * id));
    return (*Non - *Noff / (2 * id)) / sqrt(dNoff / (2 * id));
  // return (*Non - *Noff / (2 * id)) / sqrt(*Noff / (2 * id));
  else
    return -100;
}

// double gaus32(int kk, double r) {
//   int i;
//   double result = 0;
//   double reso2 = 0;

//   for (i = 0; i < 3; i++) {
//     reso2 = par[kk][i] * par[kk][i] * D2R * D2R;
//     result += par[kk][i + 3] / (2 * PI * reso2) * exp(-r * r / (2 * reso2));
//   }
//   return result;
// }

double bigaussian_mean(double theta, double sigma) {
  double thetamin = theta - 0.05 / sqrt(PI) > 0 ? theta - 0.05 / sqrt(PI) : 0;
  double thetamax = theta + 0.05 / sqrt(PI);
  // use +/- hori_sys_width/sqrt(PI) to make sure the event number in the
  // central bin is correct since the bin is a square and the integration is
  // within a circle.
  return (exp(-thetamin * thetamin / (2 * sigma * sigma)) -
          exp(-thetamax * thetamax / (2 * sigma * sigma))) /
         (2. * PI * (cos(thetamin * D2R) - cos(thetamax * D2R)) * R2D * R2D);
}

double bigaussian_simple(double theta, double sigma) {
  return exp(-theta * theta / (2 * sigma * sigma)) / (2 * PI * sigma * sigma);
}
