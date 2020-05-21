#include "constant.h"

const double D2R = 0.017453292519943296;  /*   pi/180  */
const double R2D = 57.295779513082322;    /*   180/pi  */
const double PI = 3.14159265358979323846; /* pi */

/********** Location of Tibet **********/

const double observer_la = 29.357656306;  /*  latitude  */
const double observer_lo = 100.138794639; /* longtitude */

const double kArea = 2000. * 2000.;
const int kNType = 6;
std::vector<double> kEnergyBin = {0.01, 0.1, 1., 1.e1, 1.e2, 1.e3};

const double kZenBinWidth = 0.5;
const double kZenRange = 60.;
const int kNZenBin = kZenRange / kZenBinWidth;

const int kNSmallbin = 4;
const int kPmtBinWidth = 10;
const int kPmtRange = 900;
const int kNPmtBin = kPmtRange / kPmtBinWidth;
